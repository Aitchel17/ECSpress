function info = tablegeneration_fwhmrelation(central_dir, out_dir, vessel_prefix, depth_thr, session_type, drop_flagged, max_pxsize)
%TABLEGENERATION_FWHMRELATION  Centralized FWHM heatmaps -> the tables the heatmap figures read.
%   Reads session and state out of centralized_paxfwhm_heatmap.mat (heatmap_paxfwhm samples
%   them), keeps the sessions the caller's filters pass, judges them, pools them by vessel
%   and fits the bend, and writes fwhmrelation.mat into out_dir: session, state,
%   vessel_pool, state_pool and vesselbend, with param and filt inside. The copy it
%   replaces is kept beside it as fwhmrelation_backup_<stamp>.mat.
%   Caller: e_tablegeneration
%
%   Which sessions enter is the caller's: e_tablegeneration sets the five session
%   filters. The quality gates -- which maps are drawn, which columns pool, which fits
%   count -- are this method's and are set in the body.
%
% IN   central_dir    1 x 1 str     where the centralized_*.mat files live
%      out_dir        1 x 1 str     the dataset folder fwhmrelation.mat is written to
%      vessel_prefix  1 x 1 str     VesselID prefix kept ("PA" = arteries)
%      depth_thr      1 x 1 um      kept shallower than this
%      session_type   1 x 1 str     session type kept, "" = every type
%      drop_flagged   1 x 1 bool    drop a folder whose name says notanalyzable
%      max_pxsize     1 x 1 um/px   coarsest pixel size kept, Inf = every size
% OUT  info           1 x 1 struct  n_session, out_path

arguments
    central_dir (1,1) string
    out_dir (1,1) string
    vessel_prefix (1,1) string
    depth_thr (1,1) double
    session_type (1,1) string
    drop_flagged (1,1) logical
    max_pxsize (1,1) double
end

%% Parameter

% the caller's session filters, kept in filt so fwhmrelation.mat records them
filt.vessel_prefix = vessel_prefix;
filt.depth_thr = depth_thr;
filt.session_type = session_type;
filt.drop_flagged = drop_flagged;
filt.max_pxsize = max_pxsize;
% quality gates, this method's own
filt.min_length = 150;         % Minimum data length (count), REM sets it, see FINDINGS.md
filt.min_bvamp = 2;            % Minimum diameter span (um)
filt.max_scatter = Inf;        % Maximum PVS scatter at fixed diameter (um), Inf = off
filt.pool_pxsize = [];         % Pixel size required for pooling (um/px), [] = every size
filt.pxsize_tol = 0.01;        % Pixel size tolerance (um/px)
filt.min_vessel = 10;          % Minimum vessels behind a diameter column (count)
filt.pool_common_vessel = true;% Pool every state on the vessels all of them have (logical)
filt.min_column_session = 5;   % Minimum columns per side, one session (count)
filt.min_column_vessel = 4;    % Minimum columns per side, one vessel (count)

% pooling settings; the sampling ones (y_series, event_half, state_name ...) come with the
% heatmaps and are added below
param.pool_step = [];          % Pooling grid spacing (um), [] = the coarsest pixel
param.column_stat = "median";  % Column statistic the vessel_pool curve uses (str)
param.stat_name = ["mode", "centroid", "median"];
param.state_pooled = ["awake", "drowsy", "nrem", "rem"];   % states that get a pooled map (str)
% filt forms a verdict over rows; param changes only the numbers in them

%% Read the heatmaps
heatmap_path = fullfile(central_dir, 'centralized_paxfwhm_heatmap.mat');
session_all = load(heatmap_path, 'session').session;
state_all = load(heatmap_path, 'state').state;
sampling_param = load(heatmap_path, 'param').param;
for name = string(fieldnames(sampling_param))'
    param.(name) = sampling_param.(name);
end
fprintf('centralized heatmaps : %d sessions\n   %s\n', height(session_all), heatmap_path);

%% Which sessions the result is about
cond.is_layer1 = session_all.NumericDepth < filt.depth_thr;
cond.vessel_id = string(session_all.VesselID);
cond.is_artery = startsWith(cond.vessel_id, filt.vessel_prefix, 'IgnoreCase', true);
cond.has_resolution = isfinite(session_all.NumericResolution);
% one session twice as coarse as the rest would set the pooling grid for everybody,
% and putting it on a finer grid instead would be inventing resolution it never had
cond.fine_enough = session_all.NumericResolution <= filt.max_pxsize;
% a whisker-stimulation recording and a sleep recording are different experiments,
% so two of them on one vessel are not two looks at the same thing
cond.right_type = filt.session_type == "" | ...
    string(session_all.SessionType) == filt.session_type;
cond.not_flagged = ~filt.drop_flagged | ...
    ~contains(string(session_all.Directory), "notanalyzable", 'IgnoreCase', true);
cond.keep_session = cond.is_artery & cond.is_layer1 & cond.has_resolution & ...
    cond.fine_enough & cond.right_type & cond.not_flagged;
session = session_all(cond.keep_session, :);

% message
fprintf('fwhmrelation.session : %d of %d rows kept\n', ...
    sum(cond.keep_session), numel(cond.keep_session));
fprintf('   PA %d | depth < %d um %d | resolution known %d\n', ...
    sum(cond.is_artery), filt.depth_thr, sum(cond.is_layer1), sum(cond.has_resolution));
fprintf('   type %s %d | not flagged %d | pixel <= %g um %d\n', filt.session_type, ...
    sum(cond.right_type), sum(cond.not_flagged), filt.max_pxsize, sum(cond.fine_enough));

%% Judge each session's map
% the sample-level slope is px/px = um/um, so it needs no conversion
session.enough_sample = session.n_sample >= filt.min_length;
session.wide_enough = session.enough_sample & session.bv_range >= filt.min_bvamp;
session.drawn = session.wide_enough & session.cond_sd <= filt.max_scatter;
% resolutions may be mixed: pool_by_vessel puts every map on one grid, so a pixel
% size is a reason to look at a map, not a reason to drop it
if isempty(filt.pool_pxsize)
    on_pxsize = true(height(session), 1);
else
    on_pxsize = abs(session.NumericResolution - filt.pool_pxsize) <= filt.pxsize_tol;
end
session.in_pool = session.wide_enough & on_pxsize;
% the map fits both sides whenever a line can be drawn at all; how many columns a
% side must have before that line is worth quoting is decided HERE, on counts the
% map returned, so lowering it is a re-filter and not a rebuild
session.bend_fitted = session.n_constricted >= filt.min_column_session & ...
    session.n_dilated >= filt.min_column_session;
result.fwhmrelation.session = session;


fprintf('   %d under %d samples | %d spanned less than %g um | %d drawn | %d vessel_pool\n', ...
    sum(~session.enough_sample), filt.min_length, ...
    sum(session.enough_sample & ~session.wide_enough), filt.min_bvamp, ...
    sum(session.drawn), sum(session.in_pool));
fprintf('   both sides reach %d columns : %d of %d\n', filt.min_column_session, ...
    sum(session.bend_fitted), height(session));
fprintf('\nmode-curve slope d(%s)/d(bv)\n', param.y_series);
report_slopes('spans enough diameter', session.mode_slope(session.wide_enough));


%% The same measurement again, one map per sleep state, for the sessions kept above
state = state_all(ismember(state_all.session_key, session.session_key), :);

% the same verdicts as the session table, on the same settings, applied to nothing
state.enough_sample = state.n_sample >= filt.min_length;
state.wide_enough = state.enough_sample & state.bv_range >= filt.min_bvamp;
state.drawn = state.wide_enough & state.cond_sd <= filt.max_scatter;
state.bend_fitted = state.n_constricted >= filt.min_column_session & ...
    state.n_dilated >= filt.min_column_session;
result.fwhmrelation.state = state;

fprintf('\nfwhmrelation.state : %d rows from %d sessions\n', height(state), ...
    numel(unique(state.session_key)));
for name_idx = 1:numel(param.state_name)
    on_state = state.state == param.state_name(name_idx);
    fprintf('   %-6s %2d sessions | samples %d..%d | %d span enough | %d fit both sides\n', ...
        param.state_name(name_idx), sum(on_state), min(state.n_sample(on_state)), ...
        max(state.n_sample(on_state)), sum(state.wide_enough & on_state), ...
        sum(state.bend_fitted & on_state));
end


%% Pool at the vessel level
in_pool = session(session.in_pool, :);
vessel_pool = pool_vessel(in_pool.heat, in_pool.vessel_key, in_pool.NumericResolution, ...
    param.pool_step, filt.min_vessel, param.column_stat, filt.min_column_vessel);
result.fwhmrelation.vessel_pool = vessel_pool;

%% Pool each arousal state at the vessel level
% The grid and the vessel list are the block above's, so a state map and the
% all-state map share an abscissa and a vessel index. The NORMALISATION is not --
% see pool_state_by_vessel for why a state map cannot use the column rule.
%   steady states only : a _trans interval shares frames with its neighbours. The
%   origin is the session median taken before the state mask, so the shift between
%   states is a measurement, not an offset
keep_pooled = state.drawn & ismember(state.vessel_key, vessel_pool.vessel_name);
state_drawn = state(keep_pooled, :);
fprintf('\nfwhmrelation.state_pool : %d of %d state rows kept\n', ...
    sum(keep_pooled), numel(keep_pooled));

% One state reaching further than another is the result. Left on their own vessel
% sets that comparison is confounded : a state carried by fewer vessels loses its
% outer columns to filt.min_vessel, so its curve is cut where the others are not.
if filt.pool_common_vessel
    common_vessel = unique(state_drawn.vessel_key(state_drawn.state == param.state_pooled(1)));
    for name_idx = 2:numel(param.state_pooled)
        in_this_state = unique(state_drawn.vessel_key(state_drawn.state == param.state_pooled(name_idx)));
        common_vessel = intersect(common_vessel, in_this_state);
    end
    is_common = ismember(state_drawn.vessel_key, common_vessel);
    fprintf('   on the %d vessels every state has, of %d : %d of %d rows\n', ...
        numel(common_vessel), numel(unique(state_drawn.vessel_key)), ...
        sum(is_common), numel(is_common));
    state_drawn = state_drawn(is_common, :);
end
state_pool = pool_state(state_drawn.heat, state_drawn.vessel_key, state_drawn.state, ...
    state_drawn.bv_sessionmedian, state_drawn.y_sessionmedian, vessel_pool.grid_x, ...
    vessel_pool.grid_y, vessel_pool.vessel_name, param.state_pooled, filt.min_vessel, param.column_stat);
result.fwhmrelation.state_pool = state_pool;

%% The bend, once per vessel, so it has a distribution and not just a value
vesselbend = fit_vesselbend(vessel_pool.vessel_map, vessel_pool.grid_x, vessel_pool.grid_y, ...
    vessel_pool.vessel_name, in_pool.vessel_key, in_pool.bv_sessionmedian, ...
    in_pool.eps_sessionmedian, param.stat_name, filt.min_column_vessel, param.column_stat);
result.fwhmrelation.vesselbend = vesselbend;
% the abscissa curve_* sits on. Beside the table rather than inside it, because it is
% one vector shared by every row and by vessel_pool
result.fwhmrelation.grid_x = vessel_pool.grid_x;
% the settings that produced the three, so a figure can label them without
% restating a threshold this file chose
result.fwhmrelation.param = param;
result.fwhmrelation.filt = filt;

fprintf('\nbend per vessel, the same fit off three different column statistics\n');
fprintf('%-10s %8s %12s %10s %10s %10s\n', 'statistic', 'vessels', 'constricted', ...
    'dilated', 'bend', 'p');
for stat = param.stat_name
    constricted_side = vesselbend.("constricted_" + stat);
    dilated_side = vesselbend.("dilated_" + stat);
    usable = isfinite(constricted_side) & isfinite(dilated_side);
    difference = dilated_side(usable) - constricted_side(usable);
    fprintf('%-10s %8d %12.3f %10.3f %+10.3f %10.4g\n', stat, sum(usable), ...
        median(constricted_side(usable)), median(dilated_side(usable)), ...
        median(difference), signrank(difference));
end

%% Write, one file per top-level field of result
% save_content is not a description, it is the CONVENTION: five readers across the
% tree open these files with load(path).save_content, tablegeneration_paxfwhm writes
% them the same way, and renaming it here would break every one of them silently --
% load() of a missing field errors nowhere near the file that changed it.
%   caution  this overwrites, and the copy on disk is the only record of which param built it
if ~isfolder(out_dir)
    mkdir(out_dir)
end
backup_stamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmm'));

for field_name = string(fieldnames(result))'
    out_path = fullfile(out_dir, field_name + '.mat');
    if isfile(out_path)
        backup_path = fullfile(out_dir, field_name + '_backup_' + backup_stamp + '.mat');
        copyfile(out_path, backup_path);
        fprintf('kept the previous %s as %s\n', field_name, backup_path);
    end
    save_content = result.(field_name);
    save(out_path, "save_content")
    fprintf('wrote %s\n', out_path);
    for inner_name = string(fieldnames(save_content))'
        inner = save_content.(inner_name);
        if istable(inner)
            fprintf('   %-14s table  %d x %d\n', inner_name, height(inner), width(inner));
        else
            fprintf('   %-14s %s\n', inner_name, class(inner));
        end
    end
end
info = struct('n_session', height(session), ...
    'out_path', fullfile(out_dir, string(fieldnames(result)) + ".mat"));
end

%% ---------------------------------------------------------------- helpers
function vessel_pool = pool_vessel(heat_cell, vessel_key, um_per_px, pool_step, min_vessel, ...
    column_stat, min_column_vessel)
%POOL_VESSEL  The sessions' maps put on one grid, normalised per column, averaged within and then across vessels.
%   IN   heat_cell          Nx1 cell     the sessions' heat structs
%        vessel_key         Nx1 str      which vessel each session is
%        um_per_px          Nx1 double   each session's scale
%        pool_step          1x1 double   grid spacing (um), [] = the coarsest pixel
%        min_vessel         1x1 double   vessels a column needs to be kept (count)
%        column_stat        1x1 str      the statistic the mode curve uses
%        min_column_vessel  1x1 double   columns a side needs to be fitted (count)
%   OUT  vessel_pool        1x1 struct   grid, collapsed map and mode curve, the per-vessel
%                                        maps under it, and the two-side fits off the curve
    % Sessions cannot be averaged as they stand. Each carries a different number of
    % samples and reaches a different distance along the diameter axis, so a plain mean
    % would let the widest-ranging sessions set the wings and the narrowest set the
    % centre, and the vessel_pool band's thickness would be a picture of which session
    % reached where. Normalising each DIAMETER COLUMN to sum to one first turns every
    % column into the distribution of PVS thickness AT that diameter, and the average
    % of those is a quantity with a meaning: where the PVS sits at this much dilation,
    % averaged over vessels.
    %   vessel and not session : a vessel with several sessions would otherwise weigh
    %   several times a singly-imaged one
    if isempty(pool_step)
        pool_step = max(um_per_px);
    end
    x_edge = 0;
    y_edge = 0;
    for k = 1:numel(heat_cell)
        x_edge = max(x_edge, max(abs(heat_cell{k}.x_baseceneters)));
        y_edge = max(y_edge, max(abs(heat_cell{k}.y_baseceneters)));
    end
    grid_x = -x_edge : pool_step : x_edge;
    grid_y = -y_edge : pool_step : y_edge;
    fprintf('\npooling on a %.3f um grid, %d x %d\n', pool_step, numel(grid_y), numel(grid_x));

    [vessel_name, ~, vessel_of] = unique(vessel_key, 'stable');
    [vessel_map, vessel_count] = pool_by_vessel(heat_cell, vessel_of, ...
        numel(vessel_name), grid_x, grid_y);

    vessel_reaches = sum(vessel_count > 0, 1);
    thin_column = vessel_reaches < min_vessel;
    pool_map = mean(vessel_map, 3, 'omitnan');
    pool_map(:, thin_column) = NaN;
    fprintf('   %d vessels from %d mice | columns kept %d of %d\n', numel(vessel_name), ...
        numel(unique(extractBefore(vessel_name, "_"))), sum(~thin_column), numel(grid_x));
    fprintf('   kept from %+.1f to %+.1f um, where at least %d vessels reach\n', ...
        min(grid_x(~thin_column)), max(grid_x(~thin_column)), min_vessel);

    % one number per column of the vessel_pool distribution, by the same statistic the
    % per-vessel fit uses, so the two tables are reading the same thing
    pool_mode_curve = column_statistic(pool_map, grid_y, column_stat);

    % The bend, read off the vessel_pool map rather than off one session. fit_bothsides
    % truncates the longer side, which is what makes the two slopes one measurement
    % made twice; the full-reach pair is carried beside it because the two sides do
    % not reach equally far and the difference between the pairs is worth seeing.
    vessel_pool = struct();
    vessel_pool.grid_x = grid_x;
    vessel_pool.grid_y = grid_y;
    vessel_pool.map = pool_map;
    vessel_pool.mode_curve = pool_mode_curve;
    % the stage before the collapse: one map per vessel, on the same grid. vessel_count
    % says how many of that vessel sessions reached each column, so a thin column can be
    % masked per vessel the way thin_column masks the collapsed map
    vessel_pool.vessel_map = vessel_map;
    vessel_pool.vessel_count = vessel_count;
    vessel_pool.vessel_reaches = vessel_reaches;
    vessel_pool.thin_column = thin_column;
    vessel_pool.step = pool_step;
    vessel_pool.vessel_name = vessel_name;
    [vessel_pool.sym_constricted, vessel_pool.sym_dilated, sym_info] = fit_bothsides(grid_x, pool_mode_curve, ...
        min_column_vessel);
    vessel_pool.sym_span = sym_info.span;
    vessel_pool.sym_bend = sym_info.bend;
    vessel_pool.n_sym_constricted = sym_info.n_low;
    vessel_pool.n_sym_dilated = sym_info.n_high;
    vessel_pool.sym_constricted_intercept = sym_info.intercept_low;
    vessel_pool.sym_dilated_intercept = sym_info.intercept_high;
    [vessel_pool.full_constricted, vessel_pool.n_full_constricted] = whole_reach(grid_x, pool_mode_curve, -1);
    [vessel_pool.full_dilated, vessel_pool.n_full_dilated] = whole_reach(grid_x, pool_mode_curve, 1);
    vessel_pool.r_sym_constricted = side_correlation(grid_x, pool_mode_curve, -1, vessel_pool.sym_span);
    vessel_pool.r_sym_dilated = side_correlation(grid_x, pool_mode_curve, 1, vessel_pool.sym_span);

    fprintf('   over each side''s whole reach\n');
    fprintf('      constricted %.3f (%d columns) | dilated %.3f (%d columns) | bend %+.3f\n', ...
        vessel_pool.full_constricted, vessel_pool.n_full_constricted, vessel_pool.full_dilated, vessel_pool.n_full_dilated, ...
        vessel_pool.full_dilated - vessel_pool.full_constricted);
    fprintf('   over the symmetric span, +/- %.2f um\n', vessel_pool.sym_span);
    fprintf('      constricted %.3f (%d columns) | dilated %.3f (%d columns) | bend %+.3f\n', ...
        vessel_pool.sym_constricted, vessel_pool.n_sym_constricted, vessel_pool.sym_dilated, vessel_pool.n_sym_dilated, ...
        vessel_pool.sym_bend);
end

function state_pool = pool_state(heat_cell, vessel_key, state_of, bv_median, y_median, ...
    grid_x, grid_y, vessel_name, state_pooled, min_vessel, column_stat)
%POOL_STATE  One map per arousal state, pooled by vessel on vessel_pool's grid and vessel index.
%   IN   heat_cell     Nx1 cell     the state rows' heat structs
%        vessel_key    Nx1 str      which vessel each row is
%        state_of      Nx1 str      which state each row is
%        bv_median     Nx1 double   each row's session median diameter (um)
%        y_median      Nx1 double   each row's session median of the y series (um)
%        grid_x        1xX double   vessel_pool's grid
%        grid_y        1xY double
%        vessel_name   Vx1 str      vessel_pool's vessel list, the index the maps sit on
%        state_pooled  1xS str      the states that get a map
%        min_vessel    1x1 double   vessels a column needs to be kept (count)
%        column_stat   1x1 str      the statistic the mode curve uses
%   OUT  state_pool    Sx? table    counts, the two state medians, map, mode curve, reach
    n_state = numel(state_pooled);
    state_pool = table(state_pooled(:), 'VariableNames', {'state'});
    state_pool.n_session = zeros(n_state, 1);
    state_pool.n_vessel = zeros(n_state, 1);
    state_pool.n_mouse = zeros(n_state, 1);
    state_pool.bv_statemedian = nan(n_state, 1);
    state_pool.y_statemedian = nan(n_state, 1);
    state_pool.map = cell(n_state, 1);
    state_pool.mode_curve = cell(n_state, 1);
    state_pool.vessel_reaches = cell(n_state, 1);
    state_pool.thin_column = cell(n_state, 1);

    for k = 1:n_state
        on_state = state_of == state_pooled(k);
        in_vessel = vessel_key(on_state);
        [~, vessel_index] = ismember(in_vessel, vessel_name);
        [state_map, state_count] = pool_state_by_vessel(heat_cell(on_state), vessel_index, ...
            numel(vessel_name), grid_x, grid_y);
        state_reaches = sum(state_count > 0, 1);
        state_thin = state_reaches < min_vessel;
        pooled_map = mean(state_map, 3, 'omitnan');
        pooled_map(:, state_thin) = NaN;

        state_pool.n_session(k) = sum(on_state);
        state_pool.n_vessel(k) = numel(unique(in_vessel));
        state_pool.n_mouse(k) = numel(unique(extractBefore(string(in_vessel), "_")));
        state_pool.bv_statemedian(k) = median(bv_median(on_state), 'omitnan');
        state_pool.y_statemedian(k) = median(y_median(on_state), 'omitnan');
        state_pool.map{k} = pooled_map;
        state_pool.mode_curve{k} = column_statistic(pooled_map, grid_y, column_stat);
        state_pool.vessel_reaches{k} = state_reaches;
        state_pool.thin_column{k} = state_thin;

        fprintf('   %-6s %2d sessions | %2d vessels, %d mice | columns kept %d of %d\n', ...
            state_pooled(k), sum(on_state), state_pool.n_vessel(k), ...
            state_pool.n_mouse(k), sum(~state_thin), numel(grid_x));
    end
end

function vesselbend = fit_vesselbend(vessel_map, grid_x, grid_y, vessel_name, session_vessel, ...
    bv_sessionmedian, eps_sessionmedian, stat_name, min_column_vessel, column_stat)
%FIT_VESSELBEND  Each vessel's own map fitted either side of zero, once per column statistic.
%   IN   vessel_map         YxXxV double  vessel_pool's per-vessel maps
%        grid_x             1xX double    their grid
%        grid_y             1xY double
%        vessel_name        Vx1 str       vessel_pool's vessel list, the third index of vessel_map
%        session_vessel     Nx1 str       the pooled sessions' vessel keys
%        bv_sessionmedian   Nx1 double    their median diameters (um)
%        eps_sessionmedian  Nx1 double    their median outer diameters (um)
%        stat_name          1xC str       the column statistics fitted
%        min_column_vessel  1x1 double    columns a side needs to be fitted (count)
%        column_stat        1x1 str       the statistic fitted decides
%   OUT  vesselbend         Vx? table     one row per vessel: medians, span, the two slopes and the
%                                         curve per statistic, and whether column_stat's pair fitted
    % The vessel_pool map gives one number and no spread. Fitting each vessel's own map on
    % its own symmetric span makes the two slopes a paired measurement within that
    % vessel, so the difference has a null of exactly zero and the vessels can be
    % counted.
    %   each vessel gets its own symmetric span, not a common one. The two columns per
    %   statistic are the constricted and dilated sides, fixed by the sign of x
    [~, vessel_of] = ismember(session_vessel, vessel_name);
    n_vessel = numel(vessel_name);
    vesselbend = table(vessel_name, 'VariableNames', {'vessel_key'});
    vesselbend.MouseID = extractBefore(vessel_name, "_");
    vesselbend.n_session = accumarray(vessel_of, 1);
    % where this vessel sits, so the bend can be read against calibre without a join
    vesselbend.bv_vesselmedian = accumarray(vessel_of, bv_sessionmedian, ...
        [], @(x) median(x, "omitnan"));
    vesselbend.eps_vesselmedian = accumarray(vessel_of, eps_sessionmedian, ...
        [], @(x) median(x, "omitnan"));
    vesselbend.span = nan(n_vessel, 1);
    % one row per vessel on the shared grid_x, so every vessel is a line on one axis.
    % Its own no-flux curve is not stored: bv_vesselmedian and eps_vesselmedian settle it,
    % and a stored copy could drift from them
    for stat = stat_name
        vesselbend.("constricted_" + stat) = nan(n_vessel, 1);
        vesselbend.("dilated_" + stat) = nan(n_vessel, 1);
        vesselbend.("curve_" + stat) = nan(n_vessel, numel(grid_x));
    end
    for v = 1:n_vessel
        for stat = stat_name
            curve = column_statistic(vessel_map(:, :, v), grid_y, stat);
            vesselbend.("curve_" + stat)(v, :) = curve;
            [low, high, side_info] = fit_bothsides(grid_x, curve, min_column_vessel);
            vesselbend.("constricted_" + stat)(v) = low;
            vesselbend.("dilated_" + stat)(v) = high;
            vesselbend.span(v) = side_info.span;   % the coverage is common to the three
        end
    end
    vesselbend.fitted = isfinite(vesselbend.("constricted_" + column_stat)) & ...
        isfinite(vesselbend.("dilated_" + column_stat));
end

function [vessel_map, vessel_count] = pool_by_vessel(heat_cell, vessel_of, n_vessel, ...
    grid_x, grid_y)
%POOL_BY_VESSEL  Each session's column-normalised heatmap on a common grid, summed
%   within vessel, then each column divided by the sessions that reached it.
%   OUT  vessel_map    MxNxV double  NaN in a column no session of that vessel reached
%        vessel_count  VxN double    how many of its sessions reached each column
    [mesh_x, mesh_y] = meshgrid(grid_x, grid_y);
    vessel_sum = zeros(numel(grid_y), numel(grid_x), n_vessel);
    vessel_count = zeros(n_vessel, numel(grid_x));
    for k = 1:numel(heat_cell)
        source_x = heat_cell{k}.x_baseceneters;
        source_y = heat_cell{k}.y_baseceneters;
        counts = heat_cell{k}.xy_counts_clean;
        counts(~isfinite(counts)) = 0;
        column_total = sum(counts, 1);
        filled = column_total > 0;
        counts(:, filled) = counts(:, filled) ./ column_total(filled);

        on_grid = interp2(source_x, source_y, counts, mesh_x, mesh_y, 'linear', 0);
        reaches = interp1(source_x, double(filled), grid_x, 'nearest', 0) > 0.5;
        v = vessel_of(k);
        vessel_sum(:, :, v) = vessel_sum(:, :, v) + on_grid;
        vessel_count(v, reaches) = vessel_count(v, reaches) + 1;
    end
    vessel_map = nan(numel(grid_y), numel(grid_x), n_vessel);
    for v = 1:n_vessel
        reached = vessel_count(v, :) > 0;
        if ~any(reached)
            continue
        end
        vessel_map(:, reached, v) = vessel_sum(:, reached, v) ./ vessel_count(v, reached);
    end
end

function [vessel_map, vessel_count] = pool_state_by_vessel(heat_cell, vessel_of, ...
    n_vessel, grid_x, grid_y)
%POOL_STATE_BY_VESSEL  Each session's heatmap on a common grid, normalised over the
%   WHOLE map, summed within vessel, then divided by that vessel's session count.
%   normalised over the whole map, not per column as pool_by_vessel does, so WHERE a
%   state sits on the diameter axis survives; divided by the session count, not per column
%   OUT  vessel_map    MxNxV double  sums to one per vessel; NaN for a vessel with
%                                    no session in this state
%        vessel_count  VxN double    how many of its sessions reached each column
    [mesh_x, mesh_y] = meshgrid(grid_x, grid_y);
    vessel_sum = zeros(numel(grid_y), numel(grid_x), n_vessel);
    vessel_count = zeros(n_vessel, numel(grid_x));
    session_count = zeros(n_vessel, 1);
    for k = 1:numel(heat_cell)
        source_x = heat_cell{k}.x_baseceneters;
        source_y = heat_cell{k}.y_baseceneters;
        counts = heat_cell{k}.xy_counts_clean;
        counts(~isfinite(counts)) = 0;
        column_total = sum(counts, 1);
        filled = column_total > 0;

        on_grid = interp2(source_x, source_y, counts, mesh_x, mesh_y, 'linear', 0);
        % normalised AFTER the interpolation, so a session that loses weight off the
        % edge of the shared grid still enters with the same total as every other
        grid_total = sum(on_grid, 'all');
        if grid_total <= 0
            continue
        end
        on_grid = on_grid / grid_total;
        reaches = interp1(source_x, double(filled), grid_x, 'nearest', 0) > 0.5;
        v = vessel_of(k);
        vessel_sum(:, :, v) = vessel_sum(:, :, v) + on_grid;
        vessel_count(v, reaches) = vessel_count(v, reaches) + 1;
        session_count(v) = session_count(v) + 1;
    end
    vessel_map = nan(numel(grid_y), numel(grid_x), n_vessel);
    for v = 1:n_vessel
        if session_count(v) == 0
            continue
        end
        vessel_map(:, :, v) = vessel_sum(:, :, v) / session_count(v);
    end
end

function value = column_statistic(map, grid_y, stat)
%COLUMN_STATISTIC  One number per column of a heatmap, the column read as a
%   distribution that sums to one. NaN where the column holds no weight.
    value = nan(1, size(map, 2));
    for j = 1:size(map, 2)
        weight = map(:, j);
        weight(~isfinite(weight) | weight < 0) = 0;
        total = sum(weight);
        if total <= 0
            continue
        end
        weight = weight / total;
        switch stat
            case "mode"
                [~, at_row] = max(weight);
                value(j) = grid_y(at_row);
            case "centroid"
                value(j) = sum(grid_y(:) .* weight);
            otherwise
                value(j) = half_mass(grid_y, weight);
        end
    end
end

function [slope, n_column] = whole_reach(x_value, y_value, side_sign)
%WHOLE_REACH  The unequal-leverage fit fit_bothsides refuses to make, kept only so
%   the symmetric pair can be read against it.
    on_side = isfinite(y_value) & sign(x_value) == side_sign;
    n_column = sum(on_side);
    slope = NaN;
    if n_column < 2
        return
    end
    coef = polyfit(x_value(on_side), y_value(on_side), 1);
    slope = coef(1);
end

function r_value = side_correlation(x_value, y_value, side_sign, span)
%SIDE_CORRELATION  How tight one side of the symmetric span is, against the column
%   statistic the line was fitted to.
    on_side = isfinite(y_value) & sign(x_value) == side_sign & abs(x_value) <= span;
    r_value = NaN;
    if sum(on_side) < 3
        return
    end
    r_value = corr(x_value(on_side)', y_value(on_side)');
end

function value = half_mass(centre, weight)
    % Where the column's cumulative weight crosses one half, read between bins so
    % the answer is not pinned to the grid the way the brightest bin is.
    cumulative = cumsum(weight(:));
    above = find(cumulative >= 0.5, 1);
    if isempty(above)
        value = NaN;
        return
    end
    if above == 1
        value = centre(1);
        return
    end
    below_mass = cumulative(above - 1);
    step_mass = cumulative(above) - below_mass;
    if step_mass <= 0
        value = centre(above);
        return
    end
    fraction = (0.5 - below_mass) / step_mass;
    value = centre(above - 1) + fraction * (centre(above) - centre(above - 1));
end

function report_slopes(label, slope_list)
    slope_list = slope_list(isfinite(slope_list));
    if isempty(slope_list)
        return
    end
    fprintf('   %-32s %3d sessions | median %.3f | IQR %.3f to %.3f | range %.3f to %.3f\n', ...
        label, numel(slope_list), median(slope_list), prctile(slope_list, 25), ...
        prctile(slope_list, 75), min(slope_list), max(slope_list));
end
