% Make tables that heatmap figures read.
%   Reads centralized_paxfwhm_state.mat, derives the two series from the four
%   Holding session, pooled and vesselbend
%     caution  the OUTPUT lands in the dataset folder.

clc, clear
addpath('g:\03_program\01_ecspress\09_dirstruct');   % where dirs_ecspath lives
dirs_ecspath;                                        % three roots, minus zz_notinuse

param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged_igkl_igkltdt';      % or '00_igkl' / '01_igkltdt'
end


% What gets in. Every one of these forms a verdict over sessions, over vessels or
% over diameter columns, so changing one changes WHICH data the tables are about.
param.min_sample = 150;            % samples  FWHM readouts a heatmap needs behind it.
                                   %          150 at 3 Hz is 50 s, which is what lets a
                                   %          REM row in at all -- its longest is 842
param.min_bv_range = 2;            % um       diameter span a session must cover
param.max_cond_sd = Inf;           % um       PVS scatter at a fixed diameter. Inf = off
param.pool_resolution = 0.285;     % um/px    a session is pooled only at this setting
param.resolution_tol = 0.01;       % um/px    how close counts as the same setting
param.min_vessel_n = 10;           % vessels  a diameter column needs before it is kept
param.min_column_session = 5;      % columns  a side of one session's mode curve needs
param.min_column_vessel = 4;       % columns  a side of one vessel's map needs

% What is measured. These change the numbers in the tables, not which rows are in
% them.
param.y_series = "eps";            % str      against bv; the other choice is totalpvs
param.pool_step = [];              % um       grid spacing, [] = the coarsest pixel
param.column_stat = "median";      % str      which of stat_name the pooled curve uses
% Two kinds of interval sit in state_idx and they do not align the same way.
%   *_trans  a FIXED 151-sample window, and the state it runs into starts at +75 in
%            98% of 610 events. t = 0 is therefore the window start plus 75
%   a bout    awake, nrem, rem, uarousal are scored and vary in length (a
%            microarousal runs 3 to 76 samples, 10 s at the median). t = 0 is the
%            bout start, and the window has to be wider than the bout to show it
param.event_half = 75;             % samples  either side of t = 0, transitions
param.event_half_bout = 180;       % samples  either side of t = 0, scored bouts.
                                   %          60 s at 3 Hz, six times the median
                                   %          microarousal
param.stat_name = ["mode", "centroid", "median"];
param.state_name = ["awake", "nrem", "rem", "uarousal", ...
    "an_trans", "na_trans", "nr_trans", "ra_trans"];
%   str      state_idx fields the state table splits on. The four _trans are FIXED
%            151-sample windows, not scored bouts, so every event of them is the same
%            width and they overlap the steady states they run between -- a transition
%            row and its neighbours are not independent

%% Read the centralized state analysis
dirs = dirs_central();
central_path = fullfile(dirs.central, 'centralized_paxfwhm_state.mat');
loaded_central = load(central_path); central = loaded_central.central; clear loaded_central;
% message
fprintf('centralized state : %d sessions\n   %s\n', height(central), central_path);

%% Which sessions the result is about
param.depth_thr = 70;              % um       L1 / L2-3 boundary
cond.is_layer1 = central.NumericDepth < param.depth_thr;
cond.vessel_id = string(central.VesselID);
cond.is_artery = contains(cond.vessel_id, "PA", 'IgnoreCase', true);
cond.has_resolution = isfinite(central.NumericResolution);
cond.keep_session = cond.is_artery & cond.is_layer1 & cond.has_resolution; % filtering stage
kept = central(cond.keep_session, :);

% message
fprintf('fwhmrelation.session : %d of %d rows kept\n', ...
    sum(cond.keep_session), numel(cond.keep_session));
fprintf('   PA %d | depth < %d um %d | resolution known %d\n', ...
    sum(cond.is_artery), param.depth_thr, sum(cond.is_layer1), sum(cond.has_resolution));
%% One row per session
% The heatmap and every number read off it, in one pass. There is no cheap-half /
% expensive-half split inside this file any more: re-labelling a panel costs a load
% of what this wrote, not a recompute.
%   note  the sample-level slope needs no unit conversion -- it is a ratio of two
%        thicknesses on the same pixel axis, so px/px equals um/um
measured = struct('n_sample', {}, 'bv_range', {}, 'cond_sd', {}, 'sample_slope', {}, ...
    'sample_r', {}, 'mode_slope', {}, 'slope_constricted', {}, 'slope_dilated', {}, ...
    'n_constricted', {}, 'n_dilated', {}, 'beta_ratio', {}, 'beta_gap', {}, 'bv_anchor', {}, ...
    'pvs_anchor', {}, 'area_rest', {}, 'area_gap_constricted', {}, ...
    'area_gap_dilated', {}, 'area_gapfrac_constricted', {}, ...
    'area_gapfrac_dilated', {}, 'bv_ref', {}, 'bv_rest', {}, 'pvs_rest', {}, 'heat', {});
for k = 1:height(kept)
    measured(k) = sample_heatmap(kept.data{k}.idx, kept.NumericResolution(k), ...
        param.y_series, []);
end
%%
session = kept(:, ["MouseID", "Date", "SessionType", "SessionID", "VesselID", ...
    "Directory", "NumericDepth", "NumericResolution"]);
session.vessel_key = string(kept.MouseID) + "_" + string(kept.VesselID);
session.n_sample = [measured.n_sample]';
session.bv_range = [measured.bv_range]';
session.cond_sd = [measured.cond_sd]';
session.sample_slope = [measured.sample_slope]';
session.sample_r = [measured.sample_r]';
session.mode_slope = [measured.mode_slope]';
session.slope_constricted = [measured.slope_constricted]';
session.slope_dilated = [measured.slope_dilated]';
session.n_constricted = [measured.n_constricted]';
session.n_dilated = [measured.n_dilated]';
session.beta_ratio = [measured.beta_ratio]';
session.beta_gap = [measured.beta_gap]';
session.bv_anchor = [measured.bv_anchor]';
session.pvs_anchor = [measured.pvs_anchor]';
session.area_rest = [measured.area_rest]';
session.area_gap_constricted = [measured.area_gap_constricted]';
session.area_gap_dilated = [measured.area_gap_dilated]';
session.area_gapfrac_constricted = [measured.area_gapfrac_constricted]';
session.area_gapfrac_dilated = [measured.area_gapfrac_dilated]';
session.bv_ref = [measured.bv_ref]';
session.bv_rest = [measured.bv_rest]';
session.pvs_rest = [measured.pvs_rest]';
session.heat = {measured.heat}';

%%

% the five verdicts, each named where it is defined and none of them applied here
session.enough_sample = session.n_sample >= param.min_sample;
session.wide_enough = session.enough_sample & session.bv_range >= param.min_bv_range;
session.drawn = session.wide_enough & session.cond_sd <= param.max_cond_sd;
session.in_pool = session.wide_enough & ...
    abs(session.NumericResolution - param.pool_resolution) <= param.resolution_tol;
% the map fits both sides whenever a line can be drawn at all; how many columns a
% side must have before that line is worth quoting is decided HERE, on counts the
% map returned, so lowering it is a re-filter and not a rebuild
session.bend_fitted = session.n_constricted >= param.min_column_session & ...
    session.n_dilated >= param.min_column_session;
result.fwhmrelation.session = session;


fprintf('   %d under %d samples | %d spanned less than %g um | %d drawn | %d pooled\n', ...
    sum(~session.enough_sample), param.min_sample, ...
    sum(session.enough_sample & ~session.wide_enough), param.min_bv_range, ...
    sum(session.drawn), sum(session.in_pool));
fprintf('   both sides reach %d columns : %d of %d\n', param.min_column_session, ...
    sum(session.bend_fitted), height(session));
fprintf('\nmode-curve slope d(%s)/d(bv)\n', param.y_series);
report_slopes('spans enough diameter', session.mode_slope(session.wide_enough));
held_back = session.wide_enough & ~session.drawn;
%% The same measurement again, one map per sleep state
% Not per BOUT. An NREM bout runs 91 frames at the median, which is 30 s -- too few
% to fill a heatmap and far too few to span the 2 um of diameter a slope needs. A
% session's bouts POOLED reach 2352 frames at the median and 71% of the whole
% session's diameter span, so the unit is the state within a session.
%   caution  every map re-origins on its OWN modal diameter, so the resting
%            diameter each state settles at is not in x_baseceneters. It is in
%            bv_rest, which is why that field exists
state_row = struct('n_sample', {}, 'bv_range', {}, 'cond_sd', {}, 'sample_slope', {}, ...
    'sample_r', {}, 'mode_slope', {}, 'slope_constricted', {}, 'slope_dilated', {}, ...
    'n_constricted', {}, 'n_dilated', {}, 'beta_ratio', {}, 'beta_gap', {}, 'bv_anchor', {}, ...
    'pvs_anchor', {}, 'area_rest', {}, 'area_gap_constricted', {}, ...
    'area_gap_dilated', {}, 'area_gapfrac_constricted', {}, ...
    'area_gapfrac_dilated', {}, 'bv_ref', {}, 'bv_rest', {}, 'pvs_rest', {}, 'heat', {});
state_of = strings(0);
state_session = zeros(0);
event_dbv = cell(0);
for k = 1:height(kept)
    scored = kept.data{k}.state_idx;
    n_sample_total = numel(kept.data{k}.t_axis);
    for state_name = param.state_name
        if ~isfield(scored, state_name) || isempty(scored.(state_name))
            continue
        end
        in_state = interval_mask(scored.(state_name), n_sample_total);
        state_row(end+1) = sample_heatmap(kept.data{k}.idx, ...
            kept.NumericResolution(k), param.y_series, in_state); %#ok<SAGROW>
        if endsWith(state_name, "_trans")
            onset_offset = param.event_half;
            half_width = param.event_half;
        else
            onset_offset = 0;
            half_width = param.event_half_bout;
        end
        event_dbv(end+1) = {event_average(kept.data{k}.idx, scored.(state_name), ...
            onset_offset, half_width, kept.NumericResolution(k))}; %#ok<SAGROW>
        state_of(end+1) = state_name; %#ok<SAGROW>
        state_session(end+1) = k; %#ok<SAGROW>
    end
end

state = session(state_session, ["MouseID", "Date", "SessionType", "SessionID", ...
    "VesselID", "Directory", "NumericDepth", "NumericResolution", "vessel_key"]);
state.state = state_of';
state.n_sample = [state_row.n_sample]';
state.bv_range = [state_row.bv_range]';
state.bv_ref = [state_row.bv_ref]';
state.bv_rest = [state_row.bv_rest]';
state.pvs_rest = [state_row.pvs_rest]';
state.event_dbv = event_dbv';
state.cond_sd = [state_row.cond_sd]';
state.sample_slope = [state_row.sample_slope]';
state.sample_r = [state_row.sample_r]';
state.mode_slope = [state_row.mode_slope]';
state.slope_constricted = [state_row.slope_constricted]';
state.slope_dilated = [state_row.slope_dilated]';
state.n_constricted = [state_row.n_constricted]';
state.n_dilated = [state_row.n_dilated]';
state.beta_ratio = [state_row.beta_ratio]';
state.beta_gap = [state_row.beta_gap]';
state.bv_anchor = [state_row.bv_anchor]';
state.pvs_anchor = [state_row.pvs_anchor]';
state.area_rest = [state_row.area_rest]';
state.area_gap_constricted = [state_row.area_gap_constricted]';
state.area_gap_dilated = [state_row.area_gap_dilated]';
state.area_gapfrac_constricted = [state_row.area_gapfrac_constricted]';
state.area_gapfrac_dilated = [state_row.area_gapfrac_dilated]';
state.heat = {state_row.heat}';

% the same verdicts as the session table, on the same settings, applied to nothing
state.enough_sample = state.n_sample >= param.min_sample;
state.wide_enough = state.enough_sample & state.bv_range >= param.min_bv_range;
state.drawn = state.wide_enough & state.cond_sd <= param.max_cond_sd;
state.bend_fitted = state.n_constricted >= param.min_column_session & ...
    state.n_dilated >= param.min_column_session;
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
% centre, and the pooled band's thickness would be a picture of which session
% reached where. Normalising each DIAMETER COLUMN to sum to one first turns every
% column into the distribution of PVS thickness AT that diameter, and the average
% of those is a quantity with a meaning: where the PVS sits at this much dilation,
% averaged over vessels.
%   why  vessel and not session : one vessel here contributes five sessions and
%        would otherwise weigh five times what a singly-imaged vessel weighs
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
thin_column = vessel_reaches < param.min_vessel_n;
pooled_map = mean(vessel_map, 3, 'omitnan');
pooled_map(:, thin_column) = NaN;
fprintf('   %d vessels from %d mice | columns kept %d of %d\n', numel(vessel_name), ...
    numel(unique(extractBefore(vessel_name, "_"))), sum(~thin_column), numel(grid_x));
fprintf('   kept from %+.1f to %+.1f um, where at least %d vessels reach\n', ...
    min(grid_x(~thin_column)), max(grid_x(~thin_column)), param.min_vessel_n);

% one number per column of the pooled distribution, by the same statistic the
% per-vessel fit below uses, so the two tables are reading the same thing
pooled_mode = column_statistic(pooled_map, grid_y, param.column_stat);

% The bend, read off the pooled map rather than off one session. fit_bothsides
% truncates the longer side, which is what makes the two slopes one measurement
% made twice; the full-reach pair is carried beside it because the two sides do
% not reach equally far and the difference between the pairs is worth seeing.
pooled = struct();
pooled.grid_x = grid_x;
pooled.grid_y = grid_y;
pooled.pooled_map = pooled_map;
pooled.pooled_mode = pooled_mode;
pooled.vessel_reaches = vessel_reaches;
pooled.thin_column = thin_column;
pooled.pool_step = pool_step;
pooled.vessel_name = vessel_name;
[pooled.sym_constricted, pooled.sym_dilated, sym_info] = fit_bothsides(grid_x, pooled_mode, ...
    param.min_column_vessel);
pooled.sym_span = sym_info.span;
pooled.sym_bend = sym_info.bend;
pooled.n_sym_constricted = sym_info.n_low;
pooled.n_sym_dilated = sym_info.n_high;
pooled.sym_constricted_intercept = sym_info.intercept_low;
pooled.sym_dilated_intercept = sym_info.intercept_high;
[pooled.full_constricted, pooled.n_full_constricted] = whole_reach(grid_x, pooled_mode, -1);
[pooled.full_dilated, pooled.n_full_dilated] = whole_reach(grid_x, pooled_mode, 1);
pooled.r_sym_constricted = side_correlation(grid_x, pooled_mode, -1, pooled.sym_span);
pooled.r_sym_dilated = side_correlation(grid_x, pooled_mode, 1, pooled.sym_span);
result.fwhmrelation.pooled = pooled;

fprintf('   over each side''s whole reach\n');
fprintf('      constricted %.3f (%d columns) | dilated %.3f (%d columns) | bend %+.3f\n', ...
    pooled.full_constricted, pooled.n_full_constricted, pooled.full_dilated, pooled.n_full_dilated, ...
    pooled.full_dilated - pooled.full_constricted);
fprintf('   over the symmetric span, +/- %.2f um\n', pooled.sym_span);
fprintf('      constricted %.3f (%d columns) | dilated %.3f (%d columns) | bend %+.3f\n', ...
    pooled.sym_constricted, pooled.n_sym_constricted, pooled.sym_dilated, pooled.n_sym_dilated, ...
    pooled.sym_bend);

%% The bend, once per vessel, so it has a distribution and not just a value
% The pooled map gives one number and no spread. Fitting each vessel's own map on
% its own symmetric span makes the two slopes a paired measurement within that
% vessel, so the difference has a null of exactly zero and the vessels can be
% counted.
%   why  each vessel gets ITS OWN span rather than a common one : a common span
%        would be set by the vessel that reached least far and would throw away
%        most of the others. Within a vessel the two sides are still the same
%        length, which is the only symmetry the comparison needs
%   note  the two columns per statistic are the CONSTRICTED and DILATED sides,
%        fixed by the sign of x inside fit_bothsides and never by magnitude
n_vessel = numel(vessel_name);
vesselbend = table(vessel_name, 'VariableNames', {'vessel_key'});
vesselbend.MouseID = extractBefore(vessel_name, "_");
vesselbend.n_session = accumarray(vessel_of, 1);
vesselbend.span = nan(n_vessel, 1);
for stat = param.stat_name
    vesselbend.("constricted_" + stat) = nan(n_vessel, 1);
    vesselbend.("dilated_" + stat) = nan(n_vessel, 1);
end
for v = 1:n_vessel
    for stat = param.stat_name
        curve = column_statistic(vessel_map(:, :, v), grid_y, stat);
        [low, high, side_info] = fit_bothsides(grid_x, curve, param.min_column_vessel);
        vesselbend.("constricted_" + stat)(v) = low;
        vesselbend.("dilated_" + stat)(v) = high;
        vesselbend.span(v) = side_info.span;   % the coverage is common to the three
    end
end
vesselbend.fitted = isfinite(vesselbend.("constricted_" + param.column_stat)) & ...
    isfinite(vesselbend.("dilated_" + param.column_stat));
result.fwhmrelation.vesselbend = vesselbend;
% the settings that produced the three, so a figure can label them without
% restating a threshold this file chose
result.fwhmrelation.param = param;

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
%   caution  this overwrites. The product is the pipeline's own, so overwriting is
%            allowed, but the copy on disk is the only record of which param built
%            it, and a run with one setting changed is indistinguishable afterwards
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
function row = sample_heatmap(pax_idx, um_per_px, y_series, sample_mask)
%SAMPLE_HEATMAP  One heatmap over a set of samples, and what is read off it.
%   A sample is one FWHM readout: one kymograph column, one set of four boundaries,
%   one (bv, pvs) pair. That is 1:1 with an acquired frame, but what is counted and
%   fitted here is the readout, and a readout can be missing where a frame is not.
%   The heatmap measures itself -- slopes, spread, the column counts a caller needs to
%   decide whether a slope is worth having. What is added here is the sample-level
%   pair, which reads every sample rather than the binned map, and the copy of the
%   map's numbers up into the row so the table can be queried without unpacking it.
%
%   IN   pax_idx     1x1 struct   the four boundary rows
%        um_per_px   1x1 double   this session's scale
%        y_series    1x1 str      which of fwhm_thickness's series goes on y
%        sample_mask 1xT logical  which samples to read. [] = all of them. This is
%                                 what makes a sleep state a call rather than a
%                                 second function: the same measurement, fewer samples
%   OUT  row         1x1 struct   the twenty-three fields, declared on the first line
    row = struct('n_sample', 0, 'bv_range', NaN, 'cond_sd', NaN, 'sample_slope', NaN, ...
        'sample_r', NaN, 'mode_slope', NaN, 'slope_constricted', NaN, ...
        'slope_dilated', NaN, 'n_constricted', 0, 'n_dilated', 0, ...
        'beta_ratio', NaN, 'beta_gap', NaN, 'bv_anchor', NaN, 'pvs_anchor', NaN, ...
        'area_rest', NaN, ...
        'area_gap_constricted', NaN, 'area_gap_dilated', NaN, ...
        'area_gapfrac_constricted', NaN, 'area_gapfrac_dilated', NaN, ...
        'bv_ref', NaN, 'bv_rest', NaN, 'pvs_rest', NaN, 'heat', []);

    thickness = fwhm_thickness(pax_idx);
    % COLUMNS from here down. The four boundary rows are stored as rows and
    % fwhm_thickness only subtracts, so it hands that shape straight through;
    % robustfit and corr both read a column and answer a matrix, not a number, on
    % a row. Forced here, at the reader, so line_fwhm keeps the shape it has
    % always returned
    bv_px = thickness.bv(:);
    pvs_px = thickness.(y_series)(:);
    % The anchor is the WHOLE session's median and it is computed before the mask,
    % so every state of this session and the session row itself land on one axis.
    % Median rather than the modal bin: the mode is one bin's worth of counts and
    % moves with the binning, and a state's own mode is a thing being measured here
    % rather than a thing to measure against.
    both_finite = isfinite(bv_px) & isfinite(pvs_px);
    x_origin = median(bv_px(both_finite)) * um_per_px;
    y_origin = median(pvs_px(both_finite)) * um_per_px;

    good = both_finite;
    if ~isempty(sample_mask)
        good = good & sample_mask(:);
    end
    row.n_sample = sum(good);
    bv_good = bv_px(good);
    pvs_good = pvs_px(good);


    [bin_counts, x_edge, y_edge] = histcounts2(bv_good, pvs_good, ...
        'BinWidth', [1 1], 'Normalization', 'percentage');
    % straight through, with the edges scaled -- everything the map returns is then
    % already micrometres. was xy2heatmap; see CLAUDE_LOG.md
    heat = heatmap_postprocessing(bin_counts, x_edge * um_per_px, y_edge * um_per_px, ...
        x_origin, y_origin);
    row.bv_range = range(heat.x_baseceneters);

    % up into the row, so the table answers without opening the map
    row.cond_sd = heat.column_std;
    row.mode_slope = heat.mode_slope;
    row.slope_constricted = heat.slope_constricted;
    row.slope_dilated = heat.slope_dilated;
    row.n_constricted = heat.n_constricted;
    row.n_dilated = heat.n_dilated;
    % where THIS row rests, on the session's shared axis. Two states of one vessel
    % have the same anchor and different rests, and the pair (bv_rest, pvs_rest) is
    % the only thing that separates them in the diameter plane
    row.bv_rest = heat.x_mode;
    row.pvs_rest = heat.y_mode;

    % Whether the annulus keeps its area, fitted in SQUARED coordinates. The
    % area-conserving family is eps^2 - bv^2 = const, so opening it to
    % eps^2 = beta_ratio*bv^2 + c makes the null a straight line and the fit an
    % ordinary one. Differentiating gives beta = beta_ratio * (bv/eps): the measured
    % slope is that constant multiple of the no-flux slope AT EVERY DIAMETER, which
    % is why one number covers a range over which bv/eps itself runs 0.42 to 0.65.
    %   why  samples and not the mode curve : the mode jumps from column to column
    %        and a derivative of it is noise. Squaring first, robustfit does the
    %        rest, and nothing has to be binned
    %   note  eps, never y_series. With totalpvs on y the same hypothesis sits one
    %         lower on that axis, and that shift belongs to whoever draws it
    eps_px = thickness.eps(:);
    eps_good = eps_px(good);
    % same reason as mode_slope in heatmap_postprocessing: a few-second state has
    % too few samples for a robust scale
    if numel(bv_good) >= 4
        area_fit = robustfit(bv_good.^2, eps_good.^2);
        row.beta_ratio = area_fit(2);
    end

    % where this row rests, and how much annulus it has there. The anchors are the
    % SESSION's, shared by every state of it so the axes line up; area_rest is this
    % row's own, from the samples the mask left
    row.bv_anchor = x_origin;
    row.pvs_anchor = y_origin;
    row.area_rest = (pi/4) * ((median(eps_good) * um_per_px)^2 - ...
        (median(bv_good) * um_per_px)^2);

    % The annulus area along the mode curve, ABSOLUTE and column by column. Stored
    % absolute and not as a departure because two states of one session would
    % otherwise each be measured against their own zero, and the difference between
    % them -- which is the thing being asked about -- would go into the subtraction.
    % Whoever draws it picks the baseline.
    bv_column = row.bv_anchor + heat.x_baseceneters;
    y_column = row.pvs_anchor + heat.modepvs;
    if y_series == "totalpvs"
        eps_column = bv_column + y_column;
    else
        eps_column = y_column;
    end
    heat.area_curve = (pi/4) * (eps_column.^2 - bv_column.^2);
    row.heat = heat;

    % How far the area sits from keeping itself, averaged over the columns THIS row
    % reached. A mean and not an endpoint difference: the mean is already divided by
    % the interval, so a row that swung further is not automatically further from
    % no flux. The reference is the curve at x = 0, not area_rest -- the curve is
    % built on the mode and area_rest on the medians, and mixing them would put a
    % constant offset into every departure.
    %
    % THE TWO SIDES ARE SEPARATE AND MUST STAY SEPARATE. dA/dbv = (pi/2)(beta-1)bv,
    % so with beta < 1 the area falls monotonically with diameter: the departure is
    % POSITIVE where the vessel is constricted and NEGATIVE where it is dilated, and
    % one mean over both sides cancels them. Measured on this set that cancellation
    % turned +3.3% and -4.0% into -1.4%, and what survived was not the physics but
    % the reach being twice as long on the dilated side. see CLAUDE_LOG.md
    %   caution  the columns are equally spaced, so these are means over DIAMETER and
    %            not over time. A diameter the vessel barely visited weighs the same
    %            as its resting one
    %
    % THE SPLIT IS AT THIS ROW'S OWN MEDIAN, not at the session zero. The axis is
    % anchored on the session so the states can be read against each other, but a
    % state does not sit at that anchor: REM rests about 4 um above it, so splitting
    % at zero left 22 of its 32 rows with samples on one side only and no departure
    % to measure. Constricted means constricted relative to where THAT state
    % normally sits, and its no-flux hypothesis is that the area it has THERE is
    % kept. Median and not the modal bin, for the reason the anchor uses one.
    bv_ref = median(bv_good) * um_per_px - x_origin;
    row.bv_ref = bv_ref;
    area_at_rest = interp1(heat.x_baseceneters, heat.area_curve, bv_ref);
    departure = heat.area_curve - area_at_rest;
    on_constricted = heat.x_baseceneters < bv_ref;
    on_dilated = heat.x_baseceneters > bv_ref;
    row.area_gap_constricted = mean(departure(on_constricted), "omitnan");
    row.area_gap_dilated = mean(departure(on_dilated), "omitnan");
    row.area_gapfrac_constricted = row.area_gap_constricted / area_at_rest;
    row.area_gapfrac_dilated = row.area_gap_dilated / area_at_rest;

    % beta_ratio again as a DIFFERENCE, beta - beta_noflux, which is the form that
    % walks straight into the area rate: dA/dbv = (pi/2)*eps*(beta - beta_noflux).
    % beta_noflux is bv/eps and moves with diameter (0.42 to 0.65 across one panel),
    % so a difference only means something with a diameter attached -- this one is
    % at THIS row own rest. eps there is recovered exactly from the area, no second
    % estimate and no stored intercept.
    bv_at_rest = x_origin + bv_ref;
    eps_at_rest = sqrt(bv_at_rest^2 + 4 * area_at_rest / pi);
    row.beta_gap = (row.beta_ratio - 1) * bv_at_rest / eps_at_rest;

    % the only pair the map cannot give: every sample, not every bin
    if numel(bv_good) >= 4
        sample_fit = robustfit(bv_good, pvs_good);
        row.sample_slope = sample_fit(2);
        row.sample_r = corr(bv_good, pvs_good);
    end
end

function trace = event_average(pax_idx, interval, onset_offset, half_width, um_per_px)
%EVENT_AVERAGE  Lumen diameter around the moment a transition lands, averaged over
%   this row's events. Each event is measured against ITS OWN pre-event level, so a
%   vessel that sits at 8 um and one that sits at 18 um stack on the same axis and
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
