%CENTRALIZE_BOUT  One recording's bout tables for every session, from the centralized
%   scorings, written as one table. Nothing is read from or written to the data tree.
%
%   The bout step (bout_sleep or bout_puffmove) centralize_paxfwhmstate and centralize_driftstate each
%   ran for themselves, done once. A sleep score wins over an air table, as
%   centralize_scoring decides; a sleep row holds the state bouts and their
%   transitions, an analog row the puff x locomotion bouts. The two state
%   centralizers read this instead of the scorings.
%
%   Run centralize_primary first. This reads what that wrote.

%% settings
% The columns the output table carries as its identity. util_sessionkey defines
% what a session IS and builds the lookup string; this is only which of those
% columns get written out beside the result, so the two are read together.
param.key_columns  = ["MouseID", "Date", "SessionType", "SessionID"];
param.rebuild      = false;         % bool  true recomputes every session

% Bump when bout_sleep or bout_puffmove changes what a time_table holds or how a bout is
% cut. The inputs are untouched by such a change, so no stamp can see it; every
% row is recomputed, computed moves, and the two state centralizers follow it.
%   2  puff_quietend_sec : the pre-puff quiet window ends 0.1 s before the onset
%   3  speed rule out, the puff itself is the go window, long chunks out, puffs off the air channel
schema_version = 3;

addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse
dirs = dirs_central();

%% the inputs
in_path = fullfile(dirs.central, "centralized_" + ["sleep_score", "analysis_analog"] + ".mat");
[in_exists, ~, ~] = util_checkdirstat(in_path);
if ~all(in_exists)
    error('centralize_bout:noScoring', ...
        'run centralize_primary first, centralized_sleep_score.mat or centralized_analysis_analog.mat is missing');
end
score_table = load(in_path(1)).central;
analog_table = load(in_path(2)).central;
fprintf('sleep_score %d rows, analysis_analog %d\n', height(score_table), height(analog_table));

key_score = util_sessionkey(score_table);
key_analog = util_sessionkey(analog_table);

% every session either scoring names, once each, in the order the two tables hold them
session_table = [score_table(:, param.key_columns); analog_table(:, param.key_columns)];
[~, first_i] = unique(util_sessionkey(session_table), 'stable');
session_table = session_table(first_i, :);
key_session = util_sessionkey(session_table);

%% what is already computed
out_path = fullfile(dirs.central, "centralized_bout.mat");
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
scan = struct('session_i', {}, 'scoring', {}, 'data', {}, 'input_stamp', {}, ...
    'computed', {}, 'reused', {});
for k = 1:height(session_table)
    [score, primary_analog, stamp] = centralize_scoring(score_table, key_score, ...
        analog_table, key_analog, key_session(k), 0);
    if isempty(stamp)
        fprintf('   no scoring %s\n', key_session(k));
        continue
    end
    hit = find(key_previous == key_session(k), 1);
    old = previous(hit, :);
    [row, failure] = centralize_boutsession(k, score, primary_analog, stamp, old);
    if failure ~= ""
        fprintf('   FAILED %s -- %s\n', key_session(k), failure);
    end
    if ~isempty(row)
        scan(end+1) = row; %#ok<SAGROW>
    end
    if mod(k, 25) == 0
        fprintf('   %d/%d scanned, %d held\n', k, height(session_table), numel(scan));
    end
end

%% out
kept = [scan.session_i];
reused = [scan.reused];
central = session_table(kept, :);
central.scoring = [scan.scoring]';
central.data = {scan.data}';
central.input_stamp = [scan.input_stamp]';
central.computed = [scan.computed]';
save(out_path, 'central', 'schema_version', '-v7.3');

listing = dir(out_path);
fprintf('\n%d rows: %d computed, %d reused\n', height(central), nnz(~reused), nnz(reused));
fprintf('%.0f MB -> %s\n', listing.bytes/1e6, out_path);
