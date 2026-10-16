function info = bout_paxfwhm(central_dir, rebuild)
%BOUT_PAXFWHM  Every session's FWHM series split by its bouts, through state_1d, as centralized_paxfwhm_state.mat.
%   Caller: d_centralize
%   Reads centralized_paxfwhm and centralized_bout; nothing is read from or written to the data
%   tree. A row is recomputed when its paxfwhm source moved or its bout row was recomputed.
%
% IN   central_dir  1 x 1 str    where the centralized_*.mat files live
%      rebuild      1 x 1 bool   true recomputes every session
% OUT  info         1 x 1 struct  n_rows, n_computed, n_reused, out_path

% The nine series the state analysis is run over, as <property>, <field>. Changing it changes
% what the downstream tables contain.
content_names = [ ...
    "thickness",    "eps"; ...
    "thickness",    "bv"; ...
    "thickness",    "totalpvs"; ...
    "thickness",    "uppvs"; ...
    "thickness",    "downpvs"; ...
    "displacement", "uppvs"; ...
    "displacement", "downpvs"; ...
    "displacement", "upbv"; ...
    "displacement", "downbv"];

% Bump when what is taken off state_1d changes, or when any of the five analysis methods
% changes what it produces. The inputs are untouched by such a change, so no stamp can see it.
%   2  the four boundary rows are carried alongside the results
%   3  the bout tables are centralized_bout's and are no longer copied in
%   4  state windows are half-open, and only the transition windows lose the end margin
%   5  band_decomposition and peak_trough are per bout (a NaN bout is dropped, a band
%      needs two periods of its low edge in the bout); powerdensity drops NaN bouts;
%      eventlist and event, which nothing on this path fills, are no longer carried;
%      fs is taken from t_axis, not from pax_fwhm.param.fs; param carries the psd setting;
%      input_stamp is split into input_source_modified and input_bout_computed; a transition window
%      cut by the record's first or last frame is dropped instead of clipped (bout_trim is gone)
schema_version = 5;

% 1. load centralized_bout and paxfwhm.
in_path = fullfile(central_dir, "centralized_" + ["paxfwhm", "bout"] + ".mat");
[in_exists, ~, ~] = util_checkdirstat(in_path);
if ~in_exists(1)
    error('bout_paxfwhm:noPax', 'run centralize_mat first, centralized_paxfwhm.mat is missing');
end
if ~in_exists(2)
    error('bout_paxfwhm:noBout', 'run centralize_bout first, centralized_bout.mat is missing');
end
pax_table = load(in_path(1)).central;
bout_table = load(in_path(2)).central;
fprintf('paxfwhm %d rows, bout %d\n', height(pax_table), height(bout_table));

% 2. Get session key for each data table
[key_pax, key_columns] = util_sessionkey(pax_table);
key_bout = util_sessionkey(bout_table);

% 3. Load previously processed paxstate table if exist
out_path = fullfile(central_dir, "centralized_paxfwhm_state.mat");
loaded_fwhmstatetable = table();
if isfile(out_path) && ~rebuild
    stored = load(out_path, 'central', 'schema_version');
    if isfield(stored, 'schema_version') && stored.schema_version == schema_version
        loaded_fwhmstatetable = stored.central;
    else
        fprintf('schema stamp does not match %d, recomputing every session\n', schema_version);
    end
    clear stored
end
% 3.1 Session key from loaded 
key_previous = strings(0, 1);
if ~isempty(loaded_fwhmstatetable)
    key_previous = util_sessionkey(loaded_fwhmstatetable);
end

% 4. Loop through session
scan = struct('pax_i', {}, 'data', {}, 'input_source_modified', {}, 'input_bout_computed', {}, 'reused', {});
for pax_sessionidx = 1:height(pax_table)
    % 4.1 find matching session key
    matched_boutidx = find(key_bout == key_pax(pax_sessionidx), 1);
    if isempty(matched_boutidx)
        fprintf('   no bout tables %s\n', key_pax(pax_sessionidx));
        continue
    end
    % 4.2 Select bout
    bout_data = bout_table.data{matched_boutidx};
    % 4.3 Get timestamp
    source_modified = pax_table.source_modified(pax_sessionidx);
    bout_computed = bout_table.computed(matched_boutidx);
    % 4.4 Select co-existing fwhm data session
    hit = find(key_previous == key_pax(pax_sessionidx), 1);
    stay_row = loaded_fwhmstatetable(hit, :);
    % 4.5 Determine to reuse or renew session bout X fwhm
    is_reusable = false;   % reuse only when neither the paxfwhm source moved nor the bout row was recomputed
    if ~isempty(stay_row)
        is_samesource = stay_row.input_source_modified == source_modified; % Is FWHM equal?
        is_samebout = stay_row.input_bout_computed == bout_computed; % Is bout equal?
        is_reusable = is_samesource && is_samebout; % If bout, fwhm unchanged reuse
    end
    % 5 Fill the struct array
    row = struct('pax_i', pax_sessionidx, 'data', [], 'input_source_modified', source_modified, ...
        'input_bout_computed', bout_computed, 'reused', is_reusable);
    % reuse the loaded results, or split the series again
    
    if is_reusable
        row.data = stay_row.data{1};
    else
        pax_fwhm = pax_table.data{pax_sessionidx};
        try
            row.data = get_fwhmstate(pax_fwhm, bout_data.time_table, bout_data.transition_window, content_names);
        catch err
            fprintf('   FAILED %s -- %s\n', key_pax(pax_sessionidx), err.message);
            if isempty(stay_row)
                continue
            end
            % a session that will not split keeps the loaded row, and its stamp
            row.data = stay_row.data{1};
            row.input_source_modified = stay_row.input_source_modified;
            row.input_bout_computed = stay_row.input_bout_computed;
            row.reused = true;
        end
    end
    scan(end+1) = row; %#ok<AGROW>
    if mod(pax_sessionidx, 10) == 0
        fprintf('   %d/%d scanned, %d held\n', pax_sessionidx, height(pax_table), numel(scan));
    end
end

%% out
kept = [scan.pax_i];
reused = [scan.reused];
central = pax_table(kept, key_columns);
extra = intersect(["VesselID", "Depth", "Resolution", "Cohort", "Directory"], ...
    string(pax_table.Properties.VariableNames));
central = [central, pax_table(kept, extra)];

% Depth reads '70um' and Resolution '0.19um' -- both transcribed off the
% acquisition info.txt, so the number sits at the front and the unit follows it.
% Both are read as NUMBERS downstream: Depth by the layer filter, Resolution by
% the scaling that turns pixels into micrometres, and neither can be done from
% the text. Parsed here because this is the table those readers open.
central.NumericDepth = util_unit2double(pax_table.Depth(kept));
central.NumericResolution = util_unit2double(pax_table.Resolution(kept));

central.data = {scan.data}';
central.input_source_modified = [scan.input_source_modified]';   % the stamps the next run reuses this row against
central.input_bout_computed = [scan.input_bout_computed]';
save(out_path, 'central', 'schema_version', '-v7.3');

listing = dir(out_path);
fprintf('\n%d rows: %d computed, %d reused\n', ...
    height(central), nnz(~reused), nnz(reused));
fprintf('%.0f MB -> %s\n', listing.bytes/1e6, out_path);
info = struct('n_rows', height(central), 'n_computed', nnz(~reused), 'n_reused', nnz(reused), 'out_path', out_path);
end

%% ---------------------------------------------------------------- helpers
function data = get_fwhmstate(pax_fwhm, time_table, transition_window, content_names)
%GET_FWHMSTATE  One session's FWHM series split by its bouts through state_1d: the results, the axes and the boundary rows.
%   Caller: bout_paxfwhm
%   IN   pax_fwhm           1x1 struct     its centralized_paxfwhm data
%        time_table         1x1 struct     its bout tables, centralized_bout's data.time_table
%        transition_window  1x1 float s    centralized_bout's data.transition_window, 0 for an analog row
%        content_names           Cx2 str        <property>, <field> the analysis runs over
%   OUT  data               1x1 struct     the five result structs, t_axis, state_idx, param, idx
 
    paxfwhm_state = state_1d(time_table, transition_window, pax_fwhm.t_axis);
    paxfwhm_state.get_state_indices();
    % 2. run each content the row holds through the five calculations
    for c = 1:size(content_names, 1)
        property = content_names(c, 1);
        field = content_names(c, 2);
        if ~isfield(pax_fwhm, property) || ~isfield(pax_fwhm.(property), field)
            continue
        end
        name = property + "_" + field;
        series = pax_fwhm.(property).(field);
        paxfwhm_state.get_summary(name, series);
        paxfwhm_state.get_powerdensity(name, series);
        paxfwhm_state.decompose_signal(name, series);
        paxfwhm_state.get_pppt_decomposition(name);
        paxfwhm_state.get_transitionsummary(name, series);
    end
    % 3. keep the results and the axes
    data = struct();
    for f = ["state_summary", "powerdensity", "band_decomposition", "transition", "peak_trough", ...
            "t_axis", "state_idx", "param"]
        data.(f) = paxfwhm_state.(f);
    end
    % 4. keep the four boundary rows the nine series are differences of: half their size, and they cannot disagree
    data.idx = struct();
    for f = ["clean_upperBVboundary", "clean_lowerBVboundary", ...
            "clean_pvsupedge_idx", "clean_pvsdownedge_idx"]
        data.idx.(f) = pax_fwhm.idx.(f);
    end
end
