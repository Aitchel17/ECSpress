%CENTRALIZE_PRIMARY  Every session's primary products, gathered into one file each.
%   paxfwhm.mat, polarcluster.mat and roilist.mat sit one copy per session under
%   the data tree, so any question asked across sessions costs one file open per
%   session and the answer arrives minutes later. This writes the same content as
%   three tables keyed the way the dirtable is keyed, so the same question is a
%   join.
%
%   WHAT IS DROPPED, and why. The bulk of paxfwhm.mat is rotatecrop and six
%   kymographs; the bulk of polarcluster.mat is the per-cluster median stacks.
%   Downstream reads none of it -- it reads the one-dimensional series. The two
%   raw kymographs are carried anyway, in single, because the intensity work is
%   the one thing that would otherwise have to reopen the originals. Anything
%   dropped is still in the session folder, untouched.
%
%   IT REFRESHES THE SESSION LIST FIRST. A session enters the pipeline only if the
%   merged sheet has a row for its Directory, so a recording added since the last
%   scan would be invisible here and nothing would say so. merge_dirtable rescans
%   both cohorts and rebuilds the sheet, which is about a second with the metadata
%   cache warm, and it does not overwrite anything corrected by hand.
%
%   READ ONLY under the data tree. Writes only under dirs.secondary_root -- the
%   centralized files, and the sheets merge_dirtable rebuilds.
%
%   THESE FILES ARE A CACHE, with one exception. paxfwhm, polarcluster and roilist
%   are regenerable, so no backup is kept: losing centralized_*.mat costs the
%   minutes it takes to read the sources again. sleep_score is NOT -- it is hand
%   scored and exists nowhere else, so its centralized copy is a second copy and
%   worth having for that reason alone at under a kilobyte per session. Score in
%   the session folder as before; nothing here is ever the file to edit.
%
%   sleep_score_result.fig is deliberately left out: one is larger than a whole
%   paxfwhm.mat because a .fig serialises every trace it draws, openfig wants a
%   path so a copy inside a table would have to be written back out to be opened,
%   and the .png beside it already answers what the figure is for.

%% settings
param.dataset     = "merged_igkl_igkltdt";
param.products    = [ ...            % Px2 str  <file stem>, <folder under the session>
    "paxfwhm",      "primary_analysis"; ...
    "polarcluster", "primary_analysis"; ...
    "roilist",      "primary_analysis"; ...
    "sleep_score",  "peripheral"; ...
    "analysis_analog", "peripheral"];
param.rebuild     = false;          % bool     true reads every source again

% Which cohorts go into the dataset. This is the only place it is decided --
% merge_dirtable takes it as an argument rather than keeping a second copy.
param.cohorts     = ["00_igkl", "01_igkltdt"];

% Bump this by hand whenever strip_pax or strip_polar changes which fields it
% takes. The sources are untouched by such a change, so no timestamp can see it,
% and a run that did not notice would leave rows built to two different shapes in
% one table. A mismatch rebuilds that product in full.
schema_version = 1;
% What is taken out of each object is not a setting -- it is a fact about which
% fields carry the measurement and which are bulk or regenerable. It lives with
% the code that does the taking, at the bottom of this file.

dirs = util_centraldirs();
if strlength(dirs.working_root) == 0
    error('centralize_primary:noTree', ...
        'this machine carries no data tree, so there are no sources to read');
end

% analysis_analog.mat holds an object of a class that lives in 02_othersignal, and
% ECSSession needs mdfExtractLoader out of 00_mdfExtractor, so both sibling repos
% have to be on the path or the file loads as a bare uint32. util_ecspath finds
% all three and says which is missing. Only LOADING needs them; what is written
% out is a plain struct.
addpath('g:\03_program\01_ecspress\functions');   % where util_ecspath lives
util_ecspath;                                     % three roots, minus zz_notinuse

%% refresh the session list, then read it
% merge_dirtable rescans both cohorts, rewrites their sheets and joins them. It is
% the only thing that can notice a recording added since the last scan, and this
% file finds a session only if the sheet has a row for its Directory.
%   note  the sheet is read back rather than taken as a table, because writetable
%        and readtable round-trip through Excel and every other reader of this
%        dataset sees the round-tripped version
sheet_info = merge_dirtable(dirs.working_root, dirs.secondary_root, ...
    param.cohorts, param.dataset);
if ~isempty(sheet_info.collision)
    fprintf('%d session key(s) name more than one live folder -- see the lines above\n', ...
        numel(sheet_info.collision));
end
dir_table = readtable(sheet_info.sheet_path, 'VariableNamingRule', 'preserve');
n_session = height(dir_table);
fprintf('dirtable %s\n   %d rows\n', sheet_info.sheet_path, n_session);

% The sheet keeps a row for folders that were later renamed -- a suffix gets added
% (_piv, _macrophage, _notanalyzable, ...) and both rows survive the rescan, and
% the old one points at nothing. isfolder is the whole test: a folder that was
% emptied rather than deleted still fails the isfile below, so walking every
% session to look for a file inside costs minutes and changes no outcome.
alive = false(n_session, 1);
for k = 1:n_session
    alive(k) = isfolder(string(dir_table.Directory(k)));
end
fprintf('   %d rows point at a folder, %d do not\n', sum(alive), sum(~alive));

live_row = find(alive);
live_table = dir_table(live_row, :);
key_value = util_sessionkey(live_table);

if ~isfolder(dirs.central)
    mkdir(dirs.central);
    fprintf('created %s\n', dirs.central);
end

%% one product at a time, so peak memory is one table and not three
for p = 1:size(param.products, 1)
    product = param.products(p, 1);
    subfolder = param.products(p, 2);
    file_name = product + ".mat";
    out_path = fullfile(dirs.central, "centralized_" + file_name);
    fprintf('\n=== %s ===\n', product);

    % What is already there. A row is reused only if the source file still has the
    % byte count and modified time it had when the row was made, so the check is
    % one dir() per session and never a load. The schema stamp catches the other
    % way a row goes stale: the sources untouched but strip_pax or strip_polar
    % changed which fields they take, which no timestamp can see.
    previous = table();
    if isfile(out_path) && ~param.rebuild
        stored = load(out_path, 'central', 'schema_version');
        if isfield(stored, 'schema_version') && stored.schema_version == schema_version
            previous = stored.central;
        else
            fprintf('   schema stamp does not match %d, rebuilding in full\n', schema_version);
        end
        clear stored
    end
    key_previous = strings(0, 1);
    if ~isempty(previous)
        key_previous = util_sessionkey(previous);
    end

    n_live = numel(live_row);
    data_cell = cell(n_live, 1);
    source_bytes = nan(n_live, 1);
    source_modified = nan(n_live, 1);
    has_data = false(n_live, 1);
    was_reused = false(n_live, 1);
    failed = strings(0, 1);
    for i = 1:n_live
        k = live_row(i);
        source_path = fullfile(string(dir_table.Directory(k)), subfolder, file_name);
        if ~isfile(source_path)
            continue
        end
        listing = dir(source_path);
        hit = find(key_previous == key_value(i), 1);
        if ~isempty(hit) && previous.source_bytes(hit) == listing.bytes ...
                && previous.source_modified(hit) == listing.datenum
            data_cell{i} = previous.data{hit};
            source_bytes(i) = previous.source_bytes(hit);
            source_modified(i) = previous.source_modified(hit);
            has_data(i) = true;
            was_reused(i) = true;
            continue
        end
        try
            % One variable in the file means that variable IS the content;
            % sleep_score.mat holds a dozen at the top level and the whole load
            % struct is what reads it expect, so taking the first would hand on
            % whichever happened to be alphabetically first
            loaded = load(source_path);
            field_list = fieldnames(loaded);
            if numel(field_list) == 1
                content = loaded.(field_list{1});
            else
                content = loaded;
            end
            % err  a .mat holding an object whose CLASS is not on the path does
            %      not fail to load. MATLAB warns and hands back the raw
            %      serialised bytes as a uint32, or an empty double when the
            %      object was nested inside another load -- either way the row
            %      goes on and reads later as data. One HQL080 roilist saved
            %      under the retired roi class came through as a uint32 blob.
            %      Every one of these five products is a struct or an object, so
            %      a bare numeric array IS the failure. see CLAUDE_LOG.md
            if isnumeric(content) || isempty(content)
                error('centralize_primary:classMissing', ...
                    'loaded as %s %s, not an object -- its class is not on the path: %s', ...
                    mat2str(size(content)), class(content), source_path);
            end
            switch product
                case "paxfwhm"
                    data_cell{i} = strip_pax(content);
                case "polarcluster"
                    data_cell{i} = strip_polar(content);
                case "analysis_analog"
                    data_cell{i} = strip_analog(content);
                otherwise
                    data_cell{i} = content;
            end
            source_bytes(i) = listing.bytes;
            source_modified(i) = listing.datenum;
            has_data(i) = true;
        catch err
            % a source that will not load is no reason to drop what was already
            % extracted from it, so the old row stands and the failure is named
            if ~isempty(hit)
                data_cell{i} = previous.data{hit};
                source_bytes(i) = previous.source_bytes(hit);
                source_modified(i) = previous.source_modified(hit);
                has_data(i) = true;
                was_reused(i) = true;
            end
            failed(end+1) = string(dir_table.Directory(k)) + " -- " + string(err.message); %#ok<SAGROW>
        end
        if mod(i, 25) == 0
            fprintf('   %d/%d scanned, %d read, %d reused\n', ...
                i, n_live, sum(has_data & ~was_reused), sum(was_reused));
        end
    end

    % the descriptive columns are taken fresh from the sheet every run; only data
    % and its stamp are ever carried over
    central = dir_table(live_row(has_data), :);
    central.data = data_cell(has_data);
    central.source_bytes = source_bytes(has_data);
    central.source_modified = source_modified(has_data);

    % the key has to separate the rows that are actually in this table
    key_out = key_value(has_data);
    if numel(unique(key_out)) < numel(key_out)
        [counted, name] = groupcounts(key_out);
        for k = find(counted > 1)'
            fprintf('   COLLISION %s\n', name(k));
            for h = live_row(has_data & key_value == name(k))'
                [~, stem, ext] = fileparts(string(dir_table.Directory(h)));
                fprintf('      %s\n', stem + ext);
            end
        end
    end

    % Writing 2.4 GB back unchanged is the whole cost of a run where nothing moved,
    % so the save is skipped when the new table matches the old one everywhere the
    % data is not: same rows, same order, same descriptive columns, nothing
    % recomputed and nothing dropped. Comparing without the data column keeps that
    % test cheap.
    dropped = key_previous(~ismember(key_previous, key_out));
    unchanged = all(was_reused(has_data)) && isempty(dropped) && ~isempty(previous) ...
        && height(previous) == height(central) ...
        && isequal(previous(:, ~strcmp(previous.Properties.VariableNames, 'data')), ...
                   central(:, ~strcmp(central.Properties.VariableNames, 'data')));
    if unchanged
        fprintf('   %d rows, all reused and identical -- left as it is\n', height(central));
    else
        save(out_path, 'central', 'schema_version', '-v7.3');
    end

    listing = dir(out_path);
    fprintf('   %d rows: %d reused, %d read from source, %d dropped, %d failed\n', ...
        height(central), sum(was_reused), sum(has_data & ~was_reused), ...
        numel(dropped), numel(failed));
    fprintf('   %.0f MB -> %s\n', listing.bytes/1e6, out_path);
    for k = 1:numel(dropped)
        fprintf('   DROPPED %s -- no live session carries this key any more\n', dropped(k));
    end
    for k = 1:numel(failed)
        fprintf('   FAILED %s\n', failed(k));
    end
    clear central data_cell previous
end

%% ---------------------------------------------------------------- helpers
function out = strip_pax(pax_fwhm)
    % The one-dimensional series, and the two RAW kymographs.
    % Left behind: rotatecrop and mask, which are bulk nobody downstream reads,
    % and four of the six kymographs. Those four are regenerable from these two
    % plus the cutoffs already inside param, and one of them must not be used at
    % all -- kgph_*_processed clips each frame at its own 95th percentile, which
    % is narrower than the PVS crest, so the crest comes out flattened and the
    % renormalisation then calls that flat top 1.0. see CLAUDE_LOG.md
    keep = ["idx", "thickness", "displacement", "t_axis", "up_thicker", "param"];
    raw_kgph = ["kgph_pvs", "kgph_lumen"];

    out = struct();
    for f = keep
        if isprop(pax_fwhm, f) || isfield(pax_fwhm, f)
            out.(f) = pax_fwhm.(f);
        end
    end
    out.kymograph = struct();
    for f = raw_kgph
        if isfield(pax_fwhm.kymograph, f)
            out.kymograph.(f) = single(pax_fwhm.kymograph.(f));
        end
    end
end

function out = strip_analog(primary_analog)
    % The resampled traces and the normalised spectrogram; nothing at the original
    % sample rate. rawdata, emg.power/signal and ball's pre-resample traces are all
    % 1.8 M samples per field and nothing downstream reads them. Of the four arrays
    % in ecogspectrum only baseline_S is taken -- it is the log ratio that is
    % actually looked at, and S, log_norm_spectrum and awake_S are the same
    % measurement in three other framings.
    %
    % Coming out as a plain struct also drops the dependency on the analysis_analog
    % classdef, which lives in 02_othersignal rather than either repo root here.
    % Fields are guarded because a sleep recording has no airtable or ball and a
    % whisker one has no sleep scoring behind its spectrogram.
    out = struct();
    if isprop(primary_analog, 'airtable')
        out.airtable = primary_analog.airtable;
    end
    out.ball = pick_fields(primary_analog, 'ball', ...
        ["ds_fps", "rs_taxis", "resampled_ball"]);
    out.ecog = pick_fields(primary_analog, 'ecog', ...
        ["ds_fps", "resampled_ecog", "resampled_taxis"]);
    out.emg = pick_fields(primary_analog, 'emg', ...
        ["ds_fps", "rs_taxis", "resampled_power", "resampled_signal"]);

    out.ecogspectrum = struct();
    if isprop(primary_analog, 'ecog') && isfield(primary_analog.ecog, 'ecogspectrum')
        spectrum = primary_analog.ecog.ecogspectrum;
        for f = ["T", "F", "params", "movingwin"]
            if isfield(spectrum, f)
                out.ecogspectrum.(f) = spectrum.(f);
            end
        end
        % baseline_S is the log ratio against an awake baseline, so it only exists
        % where there was sleep scoring to define that baseline. A whisker
        % recording has none and keeps the raw power S instead. The two are NOT
        % the same quantity and are not even the same way round -- baseline_S is
        % frequency by time, S is time by frequency -- so they keep their own
        % names and a reader has to look at which one is there.
        if isfield(spectrum, 'baseline_S')
            out.ecogspectrum.baseline_S = single(spectrum.baseline_S);
        elseif isfield(spectrum, 'S')
            out.ecogspectrum.S = single(spectrum.S);
        end
    end
end

function out = pick_fields(container, group, keep)
    out = struct();
    if ~isprop(container, group) || isempty(container.(group))
        return
    end
    for f = keep
        if isfield(container.(group), f)
            out.(f) = container.(group).(f);
        end
    end
end

function out = strip_polar(polarcluster)
    % Everything except clust_med_bv and clust_med_csf, which are the per-cluster
    % median image stacks and hold almost the whole file. cluster_bvcsf_constdil
    % is kept: it is the same kind of thing but only the constricted and dilated
    % pair, and the contour work reads it.
    keep = ["polar_profiles", "polar_theta", "cluster_summary", ...
        "cluster_num", "clusterboundary", "filtered_clusteridx", ...
        "constricted_cluster_idx", "dilated_cluster_idx", ...
        "constricted_medianimg", "dilated_medianimg", "cluster_bvcsf_constdil"];

    out = struct();
    for f = keep
        if isfield(polarcluster, f)
            out.(f) = polarcluster.(f);
        end
    end
end
