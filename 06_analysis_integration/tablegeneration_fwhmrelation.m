% Make tables that heatmap figures read.
%   Reads centralized_paxfwhm_state.mat, derives the two series from the four
%   Holding session, vessel_pool and vesselbend
%     caution  the OUTPUT lands in the dataset folder.

clc, clear
addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse

param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged_igkl_igkltdt';      % or '00_igkl' / '01_igkltdt'
end

%% Parameter

filt.depth_thr = 70;           % Cortical depth kept, shallower than this (um)
filt.session_type = "sleep";   % Session type kept (str), "" = every type
filt.drop_flagged = true;      % Drop a folder whose name says notanalyzable (logical)
filt.min_length = 150;         % Minimum data length (count), REM sets it, see FINDINGS.md
filt.min_bvamp = 2;            % Minimum diameter span (um)
filt.max_scatter = Inf;        % Maximum PVS scatter at fixed diameter (um), Inf = off
filt.max_pxsize = 0.3;         % Coarsest pixel size kept (um/px), Inf = every size
filt.pool_pxsize = [];         % Pixel size required for pooling (um/px), [] = every size
filt.pxsize_tol = 0.01;        % Pixel size tolerance (um/px)
filt.min_vessel = 10;          % Minimum vessels behind a diameter column (count)
filt.min_column_session = 5;   % Minimum columns per side, one session (count)
filt.min_column_vessel = 4;    % Minimum columns per side, one vessel (count)

param.y_series = "eps";        % Series on y, "eps" or "totalpvs" (str)
param.pool_step = [];          % Pooling grid spacing (um), [] = the coarsest pixel
param.column_stat = "median";  % Column statistic the vessel_pool curve uses (str)
param.event_half = 75;         % Half window, transitions (count)
param.event_half_bout = 180;   % Half window, scored bouts (count)
param.stat_name = ["mode", "centroid", "median"];
param.state_name = ["awake", "nrem", "rem", "uarousal", ...
    "an_trans", "na_trans", "nr_trans", "ra_trans"];
param.state_pooled = ["awake", "nrem", "rem"];   % states that get a pooled map (str)
% filt forms a verdict over rows; param changes only the numbers in them. t = 0 is the
% window start plus event_half for a _trans interval and the bout start for a bout;
% a _trans row overlaps its two neighbours. Why +event_half: see FINDINGS.md

%% Read the centralized state analysis
dirs = dirs_central();
central_path = fullfile(dirs.central, 'centralized_paxfwhm_state.mat');
loaded_central = load(central_path); central = loaded_central.central; clear loaded_central;
% message
fprintf('centralized state : %d sessions\n   %s\n', height(central), central_path);

%% Which sessions the result is about
cond.is_layer1 = central.NumericDepth < filt.depth_thr;
cond.vessel_id = string(central.VesselID);
cond.is_artery = contains(cond.vessel_id, "PA", 'IgnoreCase', true);
cond.has_resolution = isfinite(central.NumericResolution);
% one session twice as coarse as the rest would set the pooling grid for everybody,
% and putting it on a finer grid instead would be inventing resolution it never had
cond.fine_enough = central.NumericResolution <= filt.max_pxsize;
% a whisker-stimulation recording and a sleep recording are different experiments,
% so two of them on one vessel are not two looks at the same thing
cond.right_type = filt.session_type == "" | ...
    string(central.SessionType) == filt.session_type;
cond.not_flagged = ~filt.drop_flagged | ...
    ~contains(string(central.Directory), "notanalyzable", 'IgnoreCase', true);
cond.keep_session = cond.is_artery & cond.is_layer1 & cond.has_resolution & ...
    cond.fine_enough & cond.right_type & cond.not_flagged;
filtered = central(cond.keep_session, :);

% message
fprintf('fwhmrelation.session : %d of %d rows kept\n', ...
    sum(cond.keep_session), numel(cond.keep_session));
fprintf('   PA %d | depth < %d um %d | resolution known %d\n', ...
    sum(cond.is_artery), filt.depth_thr, sum(cond.is_layer1), sum(cond.has_resolution));
fprintf('   type %s %d | not flagged %d | pixel <= %g um %d\n', filt.session_type, ...
    sum(cond.right_type), sum(cond.not_flagged), filt.max_pxsize, sum(cond.fine_enough));
%% One row per session
% The heatmap and every number read off it, in one pass. The sample-level slope is
% px/px = um/um, so it needs no conversion
gathered = cell(height(filtered), 1);
for k = 1:height(filtered)
    thickness = fwhm_thickness(filtered.data{k}.idx);
    gathered{k} = sample_heatmap(thickness.bv, thickness.eps, param.y_series, ...
        filtered.NumericResolution(k), []);
end
measured = vertcat(gathered{:});

%% Make session table
session = filtered(:, ["MouseID", "Date", "SessionType", "SessionID", "VesselID", ...
    "Directory", "NumericDepth", "NumericResolution"]);
session.vessel_key = string(filtered.MouseID) + "_" + string(filtered.VesselID);
session = spread_rows(session, measured, "heat");


% filter condition
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
held_back = session.wide_enough & ~session.drawn;


%% The same measurement again, one map per sleep state
% pool the state within session
state_row = cell(0);
state_of = strings(0);
state_session = zeros(0);
event_dbv = cell(0);
for k = 1:height(filtered)
    scored = filtered.data{k}.state_idx;
    n_sample_total = numel(filtered.data{k}.t_axis);
    thickness = fwhm_thickness(filtered.data{k}.idx);
    for state_name = param.state_name
        if ~isfield(scored, state_name) || isempty(scored.(state_name))
            continue
        end
        in_state = interval_mask(scored.(state_name), n_sample_total); 
        state_row{end+1} = sample_heatmap(thickness.bv, thickness.eps, ...
            param.y_series, filtered.NumericResolution(k), in_state); %#ok<SAGROW>
        if endsWith(state_name, "_trans")
            onset_offset = param.event_half;
            half_width = param.event_half;
        else
            onset_offset = 0;
            half_width = param.event_half_bout;
        end
        event_dbv(end+1) = {event_average(filtered.data{k}.idx, scored.(state_name), ...
            onset_offset, half_width, filtered.NumericResolution(k))}; %#ok<SAGROW>
        state_of(end+1) = state_name; %#ok<SAGROW>
        state_session(end+1) = k; %#ok<SAGROW>
    end
end
state_row = vertcat(state_row{:});

state = session(state_session, ["MouseID", "Date", "SessionType", "SessionID", ...
    "VesselID", "Directory", "NumericDepth", "NumericResolution", "vessel_key"]);
state.state = state_of';
state = spread_rows(state, state_row, "heat");
state.event_dbv = event_dbv';

% the same verdicts as the session table, on the same settings, applied to nothing
state.enough_sample = state.n_sample >= filt.min_length;
state.wide_enough = state.enough_sample & state.bv_range >= filt.min_bvamp;
state.drawn = state.wide_enough & state.cond_sd <= filt.max_scatter;
state.bend_fitted = state.n_constricted >= filt.min_column_session & ...
    state.n_dilated >= filt.min_column_session;
result.fwhmrelation.state = state;

fprintf('\nfwhmrelation.state : %d rows from %d sessions\n', height(state), ...
    numel(unique(state_session)));
for name_idx = 1:numel(param.state_name)
    on_state = state.state == param.state_name(name_idx);
    fprintf('   %-6s %2d sessions | samples %d..%d | %d span enough | %d fit both sides\n', ...
        param.state_name(name_idx), sum(on_state), min(state.n_sample(on_state)), ...
        max(state.n_sample(on_state)), sum(state.wide_enough & on_state), ...
        sum(state.bend_fitted & on_state));
end


%% Pool at the vessel level
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
in_pool = session(session.in_pool, :);
pool_step = param.pool_step;
if isempty(pool_step)
    pool_step = max(in_pool.NumericResolution);
end
x_edge = 0;
y_edge = 0;
for k = 1:height(in_pool)
    x_edge = max(x_edge, max(abs(in_pool.heat{k}.x_baseceneters)));
    y_edge = max(y_edge, max(abs(in_pool.heat{k}.y_baseceneters)));
end
grid_x = -x_edge : pool_step : x_edge;
grid_y = -y_edge : pool_step : y_edge;
fprintf('\npooling on a %.3f um grid, %d x %d\n', pool_step, numel(grid_y), numel(grid_x));

[vessel_name, ~, vessel_of] = unique(in_pool.vessel_key, 'stable');
[vessel_map, vessel_count] = pool_by_vessel(in_pool.heat, vessel_of, ...
    numel(vessel_name), grid_x, grid_y);

vessel_reaches = sum(vessel_count > 0, 1);
thin_column = vessel_reaches < filt.min_vessel;
pool_map = mean(vessel_map, 3, 'omitnan');
pool_map(:, thin_column) = NaN;
fprintf('   %d vessels from %d mice | columns kept %d of %d\n', numel(vessel_name), ...
    numel(unique(extractBefore(vessel_name, "_"))), sum(~thin_column), numel(grid_x));
fprintf('   kept from %+.1f to %+.1f um, where at least %d vessels reach\n', ...
    min(grid_x(~thin_column)), max(grid_x(~thin_column)), filt.min_vessel);

% one number per column of the vessel_pool distribution, by the same statistic the
% per-vessel fit below uses, so the two tables are reading the same thing
pool_mode_curve = column_statistic(pool_map, grid_y, param.column_stat);

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
    filt.min_column_vessel);
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
result.fwhmrelation.vessel_pool = vessel_pool;

fprintf('   over each side''s whole reach\n');
fprintf('      constricted %.3f (%d columns) | dilated %.3f (%d columns) | bend %+.3f\n', ...
    vessel_pool.full_constricted, vessel_pool.n_full_constricted, vessel_pool.full_dilated, vessel_pool.n_full_dilated, ...
    vessel_pool.full_dilated - vessel_pool.full_constricted);
fprintf('   over the symmetric span, +/- %.2f um\n', vessel_pool.sym_span);
fprintf('      constricted %.3f (%d columns) | dilated %.3f (%d columns) | bend %+.3f\n', ...
    vessel_pool.sym_constricted, vessel_pool.n_sym_constricted, vessel_pool.sym_dilated, vessel_pool.n_sym_dilated, ...
    vessel_pool.sym_bend);

%% Pool each arousal state at the vessel level
% The grid and the vessel list are the block above's, so a state map and the
% all-state map share an abscissa and a vessel index. The NORMALISATION is not --
% see pool_state_by_vessel for why a state map cannot use the column rule.
%   steady states only : a _trans interval shares frames with its neighbours. The
%   origin is the session median taken before the state mask, so the shift between
%   states is a measurement, not an offset
keep_pooled = state.drawn & ismember(state.vessel_key, vessel_name);
state_drawn = state(keep_pooled, :);
fprintf('\nfwhmrelation.state_pool : %d of %d state rows kept\n', ...
    sum(keep_pooled), numel(keep_pooled));

n_state = numel(param.state_pooled);
state_pool = table(param.state_pooled(:), 'VariableNames', {'state'});
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
    on_state = state_drawn.state == param.state_pooled(k);
    in_state = state_drawn(on_state, :);
    [~, vessel_index] = ismember(in_state.vessel_key, vessel_name);
    [state_map, state_count] = pool_state_by_vessel(in_state.heat, vessel_index, ...
        numel(vessel_name), grid_x, grid_y);
    state_reaches = sum(state_count > 0, 1);
    state_thin = state_reaches < filt.min_vessel;
    pooled_map = mean(state_map, 3, 'omitnan');
    pooled_map(:, state_thin) = NaN;

    state_pool.n_session(k) = height(in_state);
    state_pool.n_vessel(k) = numel(unique(in_state.vessel_key));
    state_pool.n_mouse(k) = numel(unique(extractBefore(string(in_state.vessel_key), "_")));
    state_pool.bv_statemedian(k) = median(in_state.bv_sessionmedian, 'omitnan');
    state_pool.y_statemedian(k) = median(in_state.y_sessionmedian, 'omitnan');
    state_pool.map{k} = pooled_map;
    state_pool.mode_curve{k} = column_statistic(pooled_map, grid_y, param.column_stat);
    state_pool.vessel_reaches{k} = state_reaches;
    state_pool.thin_column{k} = state_thin;

    fprintf('   %-6s %2d sessions | %2d vessels, %d mice | columns kept %d of %d\n', ...
        param.state_pooled(k), height(in_state), state_pool.n_vessel(k), ...
        state_pool.n_mouse(k), sum(~state_thin), numel(grid_x));
end
result.fwhmrelation.state_pool = state_pool;

%% The bend, once per vessel, so it has a distribution and not just a value
% The vessel_pool map gives one number and no spread. Fitting each vessel's own map on
% its own symmetric span makes the two slopes a paired measurement within that
% vessel, so the difference has a null of exactly zero and the vessels can be
% counted.
%   each vessel gets its own symmetric span, not a common one. The two columns per
%   statistic are the constricted and dilated sides, fixed by the sign of x
n_vessel = numel(vessel_name);
vesselbend = table(vessel_name, 'VariableNames', {'vessel_key'});
vesselbend.MouseID = extractBefore(vessel_name, "_");
vesselbend.n_session = accumarray(vessel_of, 1);
% where this vessel sits, so the bend can be read against calibre without a join
vesselbend.bv_vesselmedian = accumarray(vessel_of, in_pool.bv_sessionmedian, ...
    [], @(x) median(x, "omitnan"));
vesselbend.eps_vesselmedian = accumarray(vessel_of, in_pool.eps_sessionmedian, ...
    [], @(x) median(x, "omitnan"));
vesselbend.span = nan(n_vessel, 1);
% one row per vessel on the shared grid_x, so every vessel is a line on one axis.
% Its own no-flux curve is not stored: bv_vesselmedian and eps_vesselmedian settle it,
% and a stored copy could drift from them
for stat = param.stat_name
    vesselbend.("constricted_" + stat) = nan(n_vessel, 1);
    vesselbend.("dilated_" + stat) = nan(n_vessel, 1);
    vesselbend.("curve_" + stat) = nan(n_vessel, numel(grid_x));
end
for v = 1:n_vessel
    for stat = param.stat_name
        curve = column_statistic(vessel_map(:, :, v), grid_y, stat);
        vesselbend.("curve_" + stat)(v, :) = curve;
        [low, high, side_info] = fit_bothsides(grid_x, curve, filt.min_column_vessel);
        vesselbend.("constricted_" + stat)(v) = low;
        vesselbend.("dilated_" + stat)(v) = high;
        vesselbend.span(v) = side_info.span;   % the coverage is common to the three
    end
end
vesselbend.fitted = isfinite(vesselbend.("constricted_" + param.column_stat)) & ...
    isfinite(vesselbend.("dilated_" + param.column_stat));
result.fwhmrelation.vesselbend = vesselbend;
% the abscissa curve_* sits on. Beside the table rather than inside it, because it is
% one vector shared by every row and by vessel_pool
result.fwhmrelation.grid_x = grid_x;
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
% tree open these files with load(path).save_content, tablegeneration_main writes
% them the same way, and renaming it here would break every one of them silently --
% load() of a missing field errors nowhere near the file that changed it.
%   caution  this overwrites, and the copy on disk is the only record of which param built it
dirs.out = fullfile(dirs.secondary_root, param.dataset);
if ~isfolder(dirs.out)
    mkdir(dirs.out)
end
backup_stamp = string(datetime('now', 'Format', 'yyyyMMdd_HHmm'));

for field_name = string(fieldnames(result))'
    out_path = fullfile(dirs.out, field_name + '.mat');
    if isfile(out_path)
        backup_path = fullfile(dirs.out, field_name + '_backup_' + backup_stamp + '.mat');
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

%% ---------------------------------------------------------------- helpers
function destination = spread_rows(destination, measured_row, boxed)
%SPREAD_ROWS  One table column per field of the measured row struct.
%   IN   destination   table       the rows these fields belong to
%        measured_row  Nx1 struct  one element per row of destination
%        boxed         1xB str     fields that are not scalars, kept in a cell
%   OUT  destination   table
    % fieldnames answers a COLUMN, and a for over a column takes it in one go
    for name = string(fieldnames(measured_row))'
        if ismember(name, boxed)
            destination.(name) = {measured_row.(name)}';
        else
            destination.(name) = [measured_row.(name)]';
        end
    end
end

function row = sample_heatmap(bv_px, eps_px, y_series, um_per_px, sample_mask)
%SAMPLE_HEATMAP  One heatmap over a set of samples, and what is read off it.

%   IN   bv_px       1xT double   lumen diameter, px
%        eps_px      1xT double   outer diameter, px
%        y_series    1x1 str      which series goes on y. "eps" or "totalpvs"
%        um_per_px   1x1 double   this session's scale
%        sample_mask 1xT logical  which samples to read. [] = all of them
%   OUT  row         1x1 struct   the twenty-three fields, declared on the first line
    row = struct('n_sample', 0, 'bv_range', NaN, 'cond_sd', NaN, 'sample_slope', NaN, ...
        'sample_r', NaN, 'mode_slope', NaN, 'slope_constricted', NaN, ...
        'slope_dilated', NaN, 'n_constricted', 0, 'n_dilated', 0, ...
        'beta_ratio', NaN, 'beta_gap', NaN, 'bv_sessionmedian', NaN, 'y_sessionmedian', NaN, ...
        'eps_sessionmedian', NaN, ...
        'area_baseline_sample', NaN, ...
        'area_gap_constricted', NaN, 'area_gap_dilated', NaN, ...
        'area_gapfrac_constricted', NaN, 'area_gapfrac_dilated', NaN, ...
        'bv_heatmedian', NaN, 'baseline_bv', NaN, 'baseline_y', NaN, 'heat', []);
    
    % COLUMNS from here down. robustfit and corr both read a column and answer a
    % matrix, not a number, on a row
    bv_px = bv_px(:);
    eps_px = eps_px(:);
    pvs_px = eps_px - bv_px;

    % what the histogram bins on y, and the lumen it dropped to get there
    if y_series == "totalpvs"
        y_px = pvs_px;
        lumen_removed = 1;
    else
        y_px = eps_px;
        lumen_removed = 0;
    end

    % Filter NaN
    both_finite = isfinite(bv_px) & isfinite(eps_px);
    good = both_finite;
    if ~isempty(sample_mask)
        good = good & sample_mask(:);
    end
    bv_good = bv_px(good); y_good = y_px(good); eps_good = eps_px(good);
    row.n_sample = sum(good);   
    
    [bin_counts, x_edge, y_edge] = histcounts2(bv_good, y_good, ...
        'BinWidth', [1 1], 'Normalization', 'percentage');

    % the shared zero for every state of this session: computed before the mask
    bv_sessionmedian = median(bv_px(both_finite)) * um_per_px;
    y_sessionmedian = median(y_px(both_finite)) * um_per_px;
    % y is eps only when y_series says so, and bv/eps needs eps either way
    eps_sessionmedian = median(eps_px(both_finite)) * um_per_px;

    % the map, then the three measurements, each on its own input
    heat = heatmap_postprocessing(bin_counts, x_edge * um_per_px, y_edge * um_per_px, ...
        bv_sessionmedian, y_sessionmedian);

    % the mode curve in absolute um, and the eps it implies
    bv_heatmedian = median(bv_good) * um_per_px - bv_sessionmedian;
    bv_column = bv_sessionmedian + heat.x_baseceneters;
    y_column = y_sessionmedian + heat.mode_curve;
    eps_column = y_column + lumen_removed * bv_column;

    slope = heatmap_modeslope(heat.x_baseceneters, heat.mode_curve, 2);
    area = heatmap_area(heat.x_baseceneters, bv_column, eps_column, bv_heatmedian);
    fitted = sample_fit(bv_good, y_good, eps_good);

    row.bv_sessionmedian = bv_sessionmedian;
    row.y_sessionmedian = y_sessionmedian;
    row.eps_sessionmedian = eps_sessionmedian;
    row.bv_range = range(heat.x_baseceneters);
    row.cond_sd = heat.column_std;
    row.baseline_bv = heat.x_mode;
    row.baseline_y = heat.y_mode;
    row.bv_heatmedian = bv_heatmedian;

    row.mode_slope = slope.whole;
    row.slope_constricted = slope.constricted;
    row.slope_dilated = slope.dilated;
    row.n_constricted = slope.n_constricted;
    row.n_dilated = slope.n_dilated;

    row.beta_ratio = fitted.beta_ratio;
    row.sample_slope = fitted.slope;
    row.sample_r = fitted.r;

    row.area_gap_constricted = area.gap_constricted;
    row.area_gap_dilated = area.gap_dilated;
    row.area_gapfrac_constricted = area.gapfrac_constricted;
    row.area_gapfrac_dilated = area.gapfrac_dilated;
    row.beta_gap = (fitted.beta_ratio - 1) * area.bv_baseline / area.eps_baseline;
    % from the median sample geometry, not from the curve above
    row.area_baseline_sample = (pi/4) * ((median(eps_good) * um_per_px)^2 - ...
        (median(bv_good) * um_per_px)^2);

    heat.area_curve = area.curve;
    row.heat = heat;
end

function fitted = sample_fit(bv_good, y_good, eps_good)
%SAMPLE_FIT  The fits made on every sample rather than on the map's columns.
%   Scale free, so pixels need no conversion.
%   IN   bv_good     Nx1 double  lumen diameter, the samples that survived the mask
%        y_good      Nx1 double  the series on y, same samples
%        eps_good    Nx1 double  outer diameter, same samples
%   OUT  fitted      1x1 struct  beta_ratio  1x1  slope of eps^2 against bv^2
%                                slope       1x1  slope of y against bv
%                                r           1x1  their correlation
fitted = struct('beta_ratio', NaN, 'slope', NaN, 'r', NaN);
if numel(bv_good) < 4
    return
end

area_fit = robustfit(bv_good.^2, eps_good.^2);
fitted.beta_ratio = area_fit(2);

series_fit = robustfit(bv_good, y_good);
fitted.slope = series_fit(2);
fitted.r = corr(bv_good, y_good);
end

function trace = event_average(pax_idx, interval, onset_offset, half_width, um_per_px)
%EVENT_AVERAGE  Lumen diameter around the moment a transition lands, averaged over
%   this row's events. Each event is measured against ITS OWN pre-event level, so a
%   narrow vessel and a wide one stack on the same axis and
%   the picture is of the excursion, not of the calibre.
%
%   Events whose window runs past either end of the recording are dropped, not
%   padded -- a padded tail would read as a flat excursion.
%
%   IN   pax_idx      1x1 struct   the four boundary rows
%        interval     Nx2 double   state_idx rows, frame indices
%        onset_offset 1x1 double   where t = 0 sits inside the interval. 75 for a
%                                  transition window, 0 for a scored bout, and the
%                                  difference is not a style choice; see param
%        half_width   1x1 double   samples either side of t = 0
%        um_per_px    1x1 double
%   OUT  trace       1 x 2*half_width+1 double   um from the pre-event level.
%                                                all NaN if no event survived
    thickness = fwhm_thickness(pax_idx);
    bv_um = thickness.bv(:) * um_per_px;
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
