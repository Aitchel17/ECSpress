%CENTRALIZE_STATE  State analysis for every session, computed from the centralized
%   inputs and written as one table. Nothing is read from or written to the data
%   tree.
%
%   This is what batch_state_analysis did, moved off the session folders. It was
%   the only stage in the chain that wrote into G:, and it wrote paxfwhm_state.mat
%   beside each recording -- a derived file living inside the raw tree, rebuilt in
%   full every run because nothing recorded what it had already done.
%
%   Same computation, same methods, same order. What changed is where the inputs
%   come from (centralized_paxfwhm and centralized_sleep_score / _analysis_analog
%   rather than a walk over 53 session folders) and where the answer goes.
%
%   Run centralize_primary first. This reads what that wrote.

%% settings
param.dataset      = "merged_igkl_igkltdt";
% The columns the output table carries as its identity. util_sessionkey defines
% what a session IS and builds the lookup string; this is only which of those
% columns get written out beside the result, so the two are read together.
param.key_columns  = ["MouseID", "Date", "SessionType", "SessionID"];
param.rebuild      = false;         % bool  true recomputes every session

% The nine series the state analysis is run over, as <property>, <field>. Not a
% setting: it is the list batch_state_analysis used, and changing it changes what
% the downstream tables contain.
param.contents = [ ...
    "thickness",    "eps"; ...
    "thickness",    "bv"; ...
    "thickness",    "totalpvs"; ...
    "thickness",    "uppvs"; ...
    "thickness",    "downpvs"; ...
    "displacement", "uppvs"; ...
    "displacement", "downpvs"; ...
    "displacement", "upbv"; ...
    "displacement", "downbv"];

% Bump when what is taken off state_linefwhm changes, or when any of the five
% analysis methods changes what it produces. The inputs are untouched by such a
% change, so no stamp can see it.
schema_version = 1;

dirs.secondary_root = ['E:\OneDrive - The Pennsylvania State University\' ...
    '2023ecspress\02_secondary_analysis'];
dirs.central = fullfile(dirs.secondary_root, 'centralized');

addpath('g:\03_program\01_ecspress\functions');   % where util_ecspath lives
util_ecspath;                                     % three roots, minus zz_notinuse

%% the inputs
in_path = fullfile(dirs.central, "centralized_" + ["paxfwhm", "sleep_score", "analysis_analog"] + ".mat");
[in_exists, ~, ~] = util_checkdirstat(in_path);
if ~in_exists(1)
    error('centralize_state:noPax', 'run centralize_primary first, centralized_paxfwhm.mat is missing');
end
pax_table = load(in_path(1)).central;
score_table = table();
analog_table = table();
if in_exists(2)
    score_table = load(in_path(2)).central;
end
if in_exists(3)
    analog_table = load(in_path(3)).central;
end
fprintf('paxfwhm %d rows, sleep_score %d, analysis_analog %d\n', ...
    height(pax_table), height(score_table), height(analog_table));

key_pax = util_sessionkey(pax_table);
key_score = util_sessionkey(score_table);
key_analog = util_sessionkey(analog_table);

%% what is already computed
out_path = fullfile(dirs.central, "centralized_paxfwhm_state.mat");
previous = table();
if isfile(out_path) && ~param.rebuild
    stored = load(out_path, 'central', 'schema_version');
    if isfield(stored, 'schema_version') && stored.schema_version == schema_version
        previous = stored.central;
    else
        fprintf('schema stamp does not match %d, recomputing every session\n', schema_version);
    end
    clear stored
end
key_previous = strings(0, 1);
if ~isempty(previous)
    key_previous = util_sessionkey(previous);
end

%% one session at a time
n_session = height(pax_table);
data_cell = cell(n_session, 1);
has_data = false(n_session, 1);
was_reused = false(n_session, 1);
input_stamp = nan(n_session, 1);
skipped = strings(0, 1);
failed = strings(0, 1);
for k = 1:n_session
    % the state analysis needs a scoring: a sleep score, or the analog analysis a
    % whisker recording carries instead. Neither means there is nothing to split by
    % An analog analysis is only a scoring when it carries stimulus times. Almost
    % every recording has an analysis_analog.mat, sleep ones included, and theirs
    % has an empty air table -- running awake_integration over it asks a table with
    % no columns for Duration. So the analog counts only when its air table has rows.
    hit_score = find(key_score == key_pax(k), 1);
    hit_analog = find(key_analog == key_pax(k), 1);
    if ~isempty(hit_analog)
        analog_row = analog_table.data{hit_analog};
        if ~isfield(analog_row, 'airtable') || isempty(analog_row.airtable)
            hit_analog = [];
        end
    end
    if isempty(hit_score) && isempty(hit_analog)
        skipped(end+1) = key_pax(k); %#ok<SAGROW>
        continue
    end

    % The stamp covers the sources this row is actually built from, and only those:
    % a sleep recording that also carries an analog analysis does not read it, so a
    % change there must not force a recompute.
    stamp = pax_table.source_modified(k);
    if ~isempty(hit_score)
        stamp = stamp + score_table.source_modified(hit_score);
    else
        stamp = stamp + analog_table.source_modified(hit_analog);
    end
    input_stamp(k) = stamp;

    hit_previous = find(key_previous == key_pax(k), 1);
    if ~isempty(hit_previous) && previous.input_stamp(hit_previous) == stamp
        data_cell{k} = previous.data{hit_previous};
        has_data(k) = true;
        was_reused(k) = true;
        continue
    end

    try
        % A sleep score wins. Most recordings carry an analysis_analog.mat whether
        % or not they were whisker sessions, so passing both would run
        % awake_integration over a sleep recording's empty air table. The old
        % batch path chose with an elseif for the same reason.
        score = [];
        primary_analog = [];
        if ~isempty(hit_score)
            score = score_table.data{hit_score};
        elseif ~isempty(hit_analog)
            primary_analog = analog_table.data{hit_analog};
        end
        pax_fwhm = pax_table.data{k};

        state_integrate = state_integration(score, primary_analog);
        state_integrate.trim_to_duration(pax_fwhm.t_axis(end));
        paxfwhm_state = state_linefwhm(state_integrate);
        paxfwhm_state.get_state_indices(pax_fwhm.t_axis, pax_fwhm.param.fs);

        for c = 1:size(param.contents, 1)
            property = param.contents(c, 1);
            field = param.contents(c, 2);
            if ~isfield(pax_fwhm, property) || ~isfield(pax_fwhm.(property), field)
                continue
            end
            name = property + "_" + field;
            data = pax_fwhm.(property).(field);
            paxfwhm_state.get_summary(name, data);
            paxfwhm_state.get_powerdensity(name, data);
            paxfwhm_state.decompose_signal(name, data);
            paxfwhm_state.get_pppt_decomposition(name);
            paxfwhm_state.get_transitionsummary(name, data);
        end

        data_cell{k} = strip_state(paxfwhm_state);
        has_data(k) = true;
    catch err
        if ~isempty(hit_previous)
            data_cell{k} = previous.data{hit_previous};
            has_data(k) = true;
            was_reused(k) = true;
        end
        failed(end+1) = key_pax(k) + " -- " + string(err.message); %#ok<SAGROW>
    end
    if mod(k, 10) == 0
        fprintf('   %d/%d, %d computed, %d reused\n', ...
            k, n_session, sum(has_data & ~was_reused), sum(was_reused));
    end
end

%% out
central = pax_table(has_data, param.key_columns);
extra = intersect(["VesselID", "Depth", "Resolution", "Cohort", "Directory"], ...
    string(pax_table.Properties.VariableNames));
central = [central, pax_table(has_data, extra)];

% Depth reads '70um' and Resolution '0.19um' -- both transcribed off the
% acquisition info.txt, so the number sits at the front and the unit follows it.
% Both are read as NUMBERS downstream: Depth by the layer filter, Resolution by
% the scaling that turns pixels into micrometres, and neither can be done from
% the text. Parsed here because this is the table those readers open.
depth_text = pax_table.Depth(has_data);
resolution_text = pax_table.Resolution(has_data);
central.NumericDepth = parse_leadnumber(depth_text);
central.NumericResolution = parse_leadnumber(resolution_text);

central.data = data_cell(has_data);
central.input_stamp = input_stamp(has_data);
save(out_path, 'central', 'schema_version', '-v7.3');

listing = dir(out_path);
fprintf('\n%d rows: %d computed, %d reused, %d with no scoring, %d failed\n', ...
    height(central), sum(has_data & ~was_reused), sum(was_reused), ...
    numel(skipped), numel(failed));
fprintf('%.0f MB -> %s\n', listing.bytes/1e6, out_path);
for k = 1:numel(failed)
    fprintf('   FAILED %s\n', failed(k));
end

%% ---------------------------------------------------------------- helpers
function value = parse_leadnumber(text_column)
    % The number at the FRONT of each entry, whatever follows it. Splitting on the
    % unit instead needs the unit spelled the way the file spells it, and
    % Resolution is written with a micro sign whose bytes depend on how the sheet
    % was saved -- a split on it returns the whole string on a mismatch and
    % str2double then gives NaN without saying why. see CLAUDE_LOG.md
    text_column = string(text_column);
    value = nan(numel(text_column), 1);
    for k = 1:numel(text_column)
        head = regexp(text_column(k), '^[\d.]+', 'match', 'once');
        if strlength(head) > 0
            value(k) = str2double(head);
        end
    end
end

function out = strip_state(paxfwhm_state)
    % The five result structs the downstream tables nest on, plus the two event
    % fields and the axes. sleep_obj is a handle to the state_integration that
    % built this, and carrying the object would carry the whole sleep score with
    % it; only the two pieces anything outside the class reads are kept --
    % rebuild_state_analysis takes param.transition_window off it, and the time
    % tables are what the windows mean.
    out = struct();
    for f = ["state_summary", "powerdensity", "band_decomposition", "transition", ...
            "peak_trough", "eventlist", "event", "t_axis", "state_idx", "param"]
        if isprop(paxfwhm_state, f)
            out.(f) = paxfwhm_state.(f);
        end
    end
    out.sleep_obj = struct();
    if ~isempty(paxfwhm_state.sleep_obj)
        out.sleep_obj.param = paxfwhm_state.sleep_obj.param;
        out.sleep_obj.time_table = paxfwhm_state.sleep_obj.time_table;
        out.sleep_obj.info = paxfwhm_state.sleep_obj.info;
    end
end
