function info = heatmap_paxfwhm(central_dir)
%HEATMAP_PAXFWHM  Every session's FWHM heatmap, whole and per state, as centralized_paxfwhm_heatmap.mat.
%   The file holds session, one row per session, and state, one row per session and scored
%   state: the session columns, session_key, vessel_key, the numbers read off the heatmap and
%   the map itself in heat; state adds event_dbv. param records the sampling settings.
%   Nothing is filtered here. Which sessions a table is about, and the quality gates read off
%   n_sample, bv_range, cond_sd, n_constricted and n_dilated, are tablegeneration_fwhmrelation's.
%   Caller: d_centralize
%   Reads central out of centralized_paxfwhm_state and samples every session on each run;
%   nothing is read from or written to the data tree.
%
% IN   central_dir  1 x 1 str    where the centralized_*.mat files live
% OUT  info         1 x 1 struct  n_rows, out_path

% sampling settings, bump schema_version when one changes
% why t = 0 sits event_half into a _trans window: see FINDINGS.md
param.y_series = "eps";        % Series on y, "eps" or "totalpvs" (str)
param.event_half = 75;         % Half window, transitions (count)
param.event_half_bout = 180;   % Half window, scored bouts (count)
param.state_name = ["awake", "drowsy", "nrem", "rem", "uarousal", ...
    "an_trans", "na_trans", "nr_trans", "ra_trans"];
schema_version = 1;

% 1. load central of the state file
paxfwhmstate_path = fullfile(central_dir, 'centralized_paxfwhm_state.mat');
if ~isfile(paxfwhmstate_path)
    error('heatmap_paxfwhm:noState', 'run bout_paxfwhm first; %s is missing', paxfwhmstate_path);
end
central = load(paxfwhmstate_path, 'central').central;
fprintf('centralized state : %d sessions\n   %s\n', height(central), paxfwhmstate_path);

% 2. key columns
id_columns = ["MouseID", "Date", "SessionType", "SessionID", "VesselID", "Directory", ...
    "NumericDepth", "NumericResolution", "Cohort", "session_key"];
id_table = central(:, id_columns);
id_table.vessel_key = string(central.MouseID) + "_" + string(central.VesselID);
% 2.1 key columns a state row repeats
state_columns = ["MouseID", "Date", "SessionType", "SessionID", "VesselID", "Directory", ...
    "NumericDepth", "NumericResolution", "session_key", "vessel_key"];

% 3. sample every session
session_rows = cell(height(id_table), 1);   % n x 1 cell of 1-row table
state_rows = cell(height(id_table), 1);     % n x 1 cell of table, empty for a session with no scored state
for k = 1:height(id_table)
    % 3.1 diameters
    thickness = fwhm_thickness(central.idx{k});
    um_per_px = central.NumericResolution(k);
    % 3.2 generate heatmap of the whole session
    whole = heatmap_sample(thickness.bv, thickness.eps, param.y_series, um_per_px, []);
    whole_measured = struct2table(whole);
    whole_measured.heat = num2cell(whole_measured.heat);
    session_rows{k} = [id_table(k, :), whole_measured];
    % 3.3 generate heatmap and event average of each scored state
    scored = central.state_idx{k};
    n_frame = numel(central.t_axis{k});
    in_state_heat = cell(1, 0);
    state_of = strings(1, 0);
    event_dbv = cell(1, 0);
    for state_name = param.state_name
        if ~isfield(scored, state_name) || isempty(scored.(state_name))
            continue
        end
        in_state = interval_mask(scored.(state_name), n_frame);
        in_state_heat{end+1} = heatmap_sample(thickness.bv, thickness.eps, ...
            param.y_series, um_per_px, in_state); %#ok<AGROW>
        if endsWith(state_name, "_trans")
            onset_offset = param.event_half;
            half_width = param.event_half;
        else
            onset_offset = 0;
            half_width = param.event_half_bout;
        end
        event_dbv{end+1} = event_average(thickness.bv, scored.(state_name), ...
            onset_offset, half_width, um_per_px); %#ok<AGROW>
        state_of(end+1) = state_name; %#ok<AGROW>
    end
    if isempty(state_of)
        continue
    end
    % 3.4 state rows
    state_key = repmat(id_table(k, state_columns), numel(state_of), 1);
    state_key.state = state_of';
    state_measured = struct2table(vertcat(in_state_heat{:}));
    state_measured.heat = num2cell(state_measured.heat);
    fragment = [state_key, state_measured];
    fragment.event_dbv = event_dbv';
    state_rows{k} = fragment;
end

% 4. stack in central's order
session = vertcat(session_rows{:});
has_state = ~cellfun(@isempty, state_rows);
state = vertcat(state_rows{has_state});

% 5. save
out_path = fullfile(central_dir, 'centralized_paxfwhm_heatmap.mat');
out_vars = struct('session', session, 'state', state, 'param', param, 'schema_version', schema_version);
save(out_path, '-struct', 'out_vars', '-v7.3');
% 5.1 summary
listing = dir(out_path);
fprintf('\n%d sessions | %d state rows\n', height(session), height(state));
fprintf('%.0f MB -> %s\n', listing.bytes/1e6, out_path);
info = struct('n_rows', height(session), 'out_path', out_path);
end

%% ---------------------------------------------------------------- helpers
function trace = event_average(bv_px, interval, onset_offset, half_width, um_per_px)
%EVENT_AVERAGE  Lumen diameter around the moment a transition lands, averaged over the events.
%   Each event is measured from its own pre-event mean; an event whose window runs past
%   either end of the recording is dropped, not padded.
%   IN   bv_px        1xT double   lumen diameter, px
%        interval     Nx2 double   state_idx rows, frame indices
%        onset_offset 1x1 double   where t = 0 sits inside the interval. event_half for a
%                                  transition window, 0 for a scored bout, and the
%                                  difference is not a style choice; see param
%        half_width   1x1 double   samples either side of t = 0
%        um_per_px    1x1 double
%   OUT  trace        1 x 2*half_width+1 double   um from the pre-event level.
%                                                 all NaN if no event survived
    bv_um = bv_px(:) * um_per_px;
    span = -half_width : half_width;
    gathered = nan(size(interval, 1), numel(span));
    for e = 1:size(interval, 1)
        onset = interval(e, 1) + onset_offset;
        take = onset + span;
        if take(1) < 1 || take(end) > numel(bv_um)
            continue
        end
        one = bv_um(take);
        gathered(e, :) = one - mean(one(span < 0), "omitnan");
    end
    trace = mean(gathered, 1, "omitnan");
end

function frame_mask = interval_mask(interval, n_frame)
%INTERVAL_MASK  state_idx's [start end] rows as one logical over frames.
%   Frames, not samples, and it is the one place the difference shows: sleep
%   scoring is written in t_axis indices. One frame yields one FWHM sample, so
%   the mask this returns lines up with the series the caller reads.
%   IN   interval    Nx2 double   frame indices, inclusive at both ends
%        n_frame     1x1 double   how long the session is
%   OUT  frame_mask  1xT logical
    frame_mask = false(1, n_frame);
    for row_idx = 1:size(interval, 1)
        first_frame = max(1, interval(row_idx, 1));
        last_frame = min(n_frame, interval(row_idx, 2));
        frame_mask(first_frame:last_frame) = true;
    end
end
