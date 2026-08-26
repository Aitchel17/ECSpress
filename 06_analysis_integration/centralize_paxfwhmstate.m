%CENTRALIZE_PAXFWHMSTATE  State analysis for every session, computed from the
%   centralized inputs and written as one table. Nothing is read from or written
%   to the data tree.
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
%   2  the four boundary rows are carried alongside the results
schema_version = 2;

dirs = dirs_central();

addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse

%% the inputs
in_path = fullfile(dirs.central, "centralized_" + ["paxfwhm", "sleep_score", "analysis_analog"] + ".mat");
[in_exists, ~, ~] = util_checkdirstat(in_path);
if ~in_exists(1)
    error('centralize_paxfwhmstate:noPax', 'run centralize_primary first, centralized_paxfwhm.mat is missing');
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

%% one row per session that has a scoring to split by
scan = struct('pax_i', {}, 'data', {}, 'input_stamp', {}, 'reused', {});
for k = 1:height(pax_table)
    [score, primary_analog, stamp] = centralize_scoring(score_table, key_score, ...
        analog_table, key_analog, key_pax(k), pax_table.source_modified(k));
    if isempty(stamp)
        fprintf('   no scoring %s\n', key_pax(k));
        continue
    end
    hit = find(key_previous == key_pax(k), 1);
    old = previous(hit, :);
    [row, failure] = centralize_statesession(k, pax_table.data{k}, score, ...
        primary_analog, param.contents, stamp, old);
    if failure ~= ""
        fprintf('   FAILED %s -- %s\n', key_pax(k), failure);
    end
    if ~isempty(row)
        scan(end+1) = row; %#ok<SAGROW>
    end
    if mod(k, 10) == 0
        fprintf('   %d/%d scanned, %d held\n', k, height(pax_table), numel(scan));
    end
end

%% out
kept = [scan.pax_i];
reused = [scan.reused];
central = pax_table(kept, param.key_columns);
extra = intersect(["VesselID", "Depth", "Resolution", "Cohort", "Directory"], ...
    string(pax_table.Properties.VariableNames));
central = [central, pax_table(kept, extra)];

% Depth reads '70um' and Resolution '0.19um' -- both transcribed off the
% acquisition info.txt, so the number sits at the front and the unit follows it.
% Both are read as NUMBERS downstream: Depth by the layer filter, Resolution by
% the scaling that turns pixels into micrometres, and neither can be done from
% the text. Parsed here because this is the table those readers open.
central.NumericDepth = parse_leadnumber(pax_table.Depth(kept));
central.NumericResolution = parse_leadnumber(pax_table.Resolution(kept));

central.data = {scan.data}';
central.input_stamp = [scan.input_stamp]';
save(out_path, 'central', 'schema_version', '-v7.3');

listing = dir(out_path);
fprintf('\n%d rows: %d computed, %d reused\n', ...
    height(central), nnz(~reused), nnz(reused));
fprintf('%.0f MB -> %s\n', listing.bytes/1e6, out_path);

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
