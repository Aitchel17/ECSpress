%FIG_FWHM_AREAGAP  How far the annulus area sits from keeping itself, and the spread.
%   One block, one figure. Nothing is tiled and nothing is saved unless you ask, so
%   a block can be run on its own while looking at the one before it.
%
%   The unit of replication is the VESSEL. One vessel here contributes up to five
%   sessions, so sessions are averaged within vessel first and every n below counts
%   vessels. see tablegeneration_fwhmrelation
clc, clear
%% 1. Read

addpath('g:\03_program\01_ecspress\09_dirstruct');
dirs_ecspath;
clee = color_lee();

param.dataset = getenv('ECSPRESS_DATASET');
if isempty(param.dataset)
    param.dataset = 'merged_igkl_igkltdt';
end
% Name and label side by side, matched by POSITION. The table stores an
% identifier and an axis names the thing; a subset takes the same indices out of
% both, which is what param.*_pick below is for
param.state_name  = ["all",           "awake", "nrem", "rem", "uarousal", ...
                     "an_trans",      "na_trans",      "nr_trans",    "ra_trans"];
param.state_label = ["whole session", "awake", "NREM", "REM", "microarousal", ...
                     "awake to NREM", "NREM to awake", "NREM to REM", "REM to awake"];
param.scatter_pick = [2 3 4];      % the states that hold still, for the scatter
param.paired_pick  = [3 4 5 6 7 8 9];   % everything the reference is compared to
% the scatter puts one point per row, so it only takes the states that hold still.
% A transition window overlaps the two states it runs between and would draw a
% third point on top of them
param.scatter_state = param.state_name(param.scatter_pick);
param.scatter_label = param.state_label(param.scatter_pick);
param.fontsize = 11;

dirs = dirs_central();
dirs.save_dir = fullfile(dirs.secondary_root, param.dataset);
relation = load(fullfile(dirs.save_dir, 'fwhmrelation.mat')).save_content;

% the session table is the whole recording, so it is the "all" state
session = relation.session;
session.state = repmat("all", height(session), 1);
shared = intersect(session.Properties.VariableNames, relation.state.Properties.VariableNames);
rows = [session(:, shared); relation.state(:, shared)];
rows = rows(rows.wide_enough, :);
fprintf('%d rows over %d vessels\n', height(rows), numel(unique(rows.vessel_key)));

%% 2. Where the vessels sit
% bv against eps, one point per session. The iso-area curves are what makes the
% spread readable: two vessels on one curve hold the same annulus at rest however
% different their diameters are.
fig_scatter = figure('Color', 'w', 'Position', [80 80 620 580], 'Name', 'fig_fwhm_areagap_scatter');
ax = axes(fig_scatter);
hold(ax, 'on')

bv_line = linspace(4, 22, 200);
for iso_area = [100 200 400 800]
    eps_iso = sqrt(bv_line.^2 + (4/pi) * iso_area);
    plot(ax, bv_line, eps_iso, '-', 'Color', [0.85 0.85 0.85], 'LineWidth', 1);
    text(ax, bv_line(end), eps_iso(end), sprintf(' %d um^2', iso_area), ...
        'Color', [0.6 0.6 0.6], 'FontSize', param.fontsize - 3);
end
plot(ax, bv_line, bv_line, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5);

% one point per VESSEL. A vessel here carries up to five sessions and drawing each
% of them would weigh a repeatedly imaged vessel five times, which is the same
% reason every n below counts vessels
[vessel, ~, vessel_of] = unique(rows.vessel_key);
n_vessel = numel(vessel);
bv_vessel = column_mean(rows, vessel_of, n_vessel, "bv_sessionmedian", "all");
y_vessel = column_mean(rows, vessel_of, n_vessel, "y_sessionmedian", "all");
gap_vessel = column_mean(rows, vessel_of, n_vessel, "area_gapfrac_dilated", "all");
scatter(ax, bv_vessel, y_vessel, 46, gap_vessel, ...
    'filled', 'MarkerEdgeColor', [0.25 0.25 0.25]);
colormap(ax, clee.gradient.inferno)
bar_handle = colorbar(ax);
bar_handle.Label.String = 'area gap on the dilated side, fraction of rest';
xlabel(ax, 'lumen diameter at rest, um')
ylabel(ax, 'outer PVS diameter at rest, um')
set(ax, 'FontSize', param.fontsize, 'Box', 'on')
axis(ax, 'equal')
xlim(ax, [4 22])
title(ax, sprintf('%d vessels | grey = equal annulus area | diagonal = no PVS', ...
    sum(isfinite(bv_vessel))))

%% 2b. Where each STATE rests, and the line from one to the other
% The anchors are session-level, so plotting state rows at (bv_sessionmedian, y_sessionmedian)
% would stack every state of a session on one point. What separates them is where
% each one SETTLES on that shared axis. Note the two are not the same estimator --
% bv_heatmedian is a median and baseline_y is a modal ridge; see CLAUDE_LOG.md
%   the grey line joins the states of one vessel-session, so its direction is the
%   move from one state to the other and its length is how far
fig_state = figure('Color', 'w', 'Position', [100 100 620 580], 'Name', 'fig_fwhm_areagap_state');
ax = axes(fig_state);
hold(ax, 'on')

bv_line = linspace(4, 22, 200);
for iso_area = [100 200 400 800]
    plot(ax, bv_line, sqrt(bv_line.^2 + (4/pi) * iso_area), '-', ...
        'Color', [0.88 0.88 0.88], 'LineWidth', 1);
end

% from param, not a second list here. This block carried its own ["nrem","awake"]
% and so left REM out of a figure it belongs in
state_shown = param.scatter_state;
state_colour = cell2mat(arrayfun(@(f) state_rgb(clee, f), state_shown', ...
    'UniformOutput', false));

% absolute_bv and absolute_y are sums of two columns, and a mean of a sum is the
% sum of the means, so each is averaged within vessel one column at a time
absolute_bv = nan(n_vessel, numel(state_shown));
absolute_y = nan(n_vessel, numel(state_shown));
for s = 1:numel(state_shown)
    absolute_bv(:, s) = column_mean(rows, vessel_of, n_vessel, "bv_sessionmedian", state_shown(s)) + ...
        column_mean(rows, vessel_of, n_vessel, "bv_heatmedian", state_shown(s));
    absolute_y(:, s) = column_mean(rows, vessel_of, n_vessel, "y_sessionmedian", state_shown(s)) + ...
        column_mean(rows, vessel_of, n_vessel, "baseline_y", state_shown(s));
end

% one line per vessel, drawn first so the points sit on top of it
for v = 1:n_vessel
    on_v = isfinite(absolute_bv(v, :)) & isfinite(absolute_y(v, :));
    if sum(on_v) < 2
        continue
    end
    plot(ax, absolute_bv(v, on_v), absolute_y(v, on_v), '-', 'Color', [0.7 0.7 0.7], ...
        'LineWidth', 0.8);
end

handle_state = gobjects(1, numel(state_shown));
for s = 1:numel(state_shown)
    handle_state(s) = scatter(ax, absolute_bv(:, s), absolute_y(:, s), 44, ...
        state_colour(s, :), 'filled', 'MarkerEdgeColor', [0.25 0.25 0.25]);
end
legend(ax, handle_state, param.scatter_label, 'Location', 'northwest', 'Box', 'off')
xlabel(ax, 'lumen diameter where that state settles, um')
ylabel(ax, 'outer PVS diameter there, um')
set(ax, 'FontSize', param.fontsize, 'Box', 'on')
axis(ax, 'equal')
paired = sum(sum(isfinite(absolute_bv), 2) >= 2);
title(ax, sprintf('%d vessels with two or more states | grey = equal annulus area', ...
    paired))

% how far the states sit apart, on the vessels that carry them
fprintf('%-14s %10s %10s\n', 'state', 'baseline bv', 'baseline pvs');
for s = 1:numel(state_shown)
    fprintf('%-14s %10.2f %10.2f\n', param.scatter_label(s), ...
        median(absolute_bv(:, s), 'omitnan'), median(absolute_y(:, s), 'omitnan'));
end

%% 3. Area gap by state, um^2
% Bar is the mean over vessels, thick whisker the 95% CI, thin whisker one SD, and
% every vessel is drawn so the CI can be read against what it came from.
%   two bars per state : the departure is positive constricted and negative dilated,
%   and one mean over both sides cancels them. see CLAUDE_LOG.md
[gap_um, label_um] = by_state(rows, ["area_gap_constricted", "area_gap_dilated"], ...
    param.state_name, param.state_label);
fig_gap = figure('Color', 'w', 'Position', [110 110 720 520], 'Name', 'fig_fwhm_areagap_um');
draw_bars(axes(fig_gap), gap_um, label_um, param.fontsize, ...
    'area gap, um^2', 'left bar constricted, right bar dilated');

%% 4. Area gap by state, as a fraction of the resting annulus
% The same quantity with vessel size divided out: the micrometre bars above carry the
% spread of area_baseline_sample inside them. For its range see FINDINGS.md
[gap_frac, label_frac] = by_state(rows, ["area_gapfrac_constricted", ...
    "area_gapfrac_dilated"], param.state_name, param.state_label);
fig_frac = figure('Color', 'w', 'Position', [140 140 720 520], 'Name', 'fig_fwhm_areagap_frac');
draw_bars(axes(fig_frac), gap_frac, label_frac, param.fontsize, ...
    'area gap, fraction of the resting annulus', 'left bar constricted, right bar dilated');

%% 5. Every state against awake, PAIRED within vessel
% The unpaired bars average sets that are not the same vessels, so the difference
% between them carries the difference between the sets. Here only vessels holding
% BOTH states of a pair are kept and each is one line, so the comparison is within
% vessel and the vessel-to-vessel spread -- the largest thing on this figure --
% drops out. The test is signrank on those within-vessel differences.
%   awake is the reference, as in tablegeneration_main. The awake bar is recomputed
%   for each pair, on that pair's own vessels
param.paired_ref = "awake";
% beta - beta_noflux and not the ratio: the difference is what appears in
% dA/dbv = (pi/2)*eps*(beta - beta_noflux), so a bar height reads straight as an
% area rate. Its null is 0, the ratio null was 1. Set this to "beta_ratio" to see
% the other form; the vessels rank the same either way (r +0.98, 30 of 30 signs)
param.beta_column = "beta_gap";
param.beta_null = 0;
param.beta_label = "beta minus the no-flux beta, at rest";
param.paired_other = param.state_name(param.paired_pick);
param.paired_other_label = param.state_label(param.paired_pick);
param.paired_ref_label = "awake";
param.paired_column = ["area_gapfrac_constricted", "area_gapfrac_dilated"];
param.paired_label = ["constricted", "dilated"];

fig_paired = figure('Color', 'w', 'Position', [170 170 820 540], 'Name', 'fig_fwhm_areagap_paired');
draw_paired(axes(fig_paired), rows, param.paired_column, param.paired_label, ...
    param.paired_ref, param.paired_ref_label, param.paired_other, ...
    param.paired_other_label, param.fontsize, ...
    'area gap, fraction of the resting annulus');

%% 5c. The same difference, with the flagged vessels drawn apart
% A vessel whose AWAKE fit lands outside 0..1 is doing something the model does not
% describe: below 0 the outer boundary moves against the lumen, above 1 the annulus
% gains area as the vessel dilates. They are drawn as their own group rather than
% dropped -- if the two groups behave the same the flag does not matter.
%   caution  the flag is on AWAKE only, so it cannot select on the thing being compared
param.flag_beta_range = [0 1];   % awake beta_ratio inside this is not flagged

awake_beta = rows(rows.state == "awake", ["vessel_key", "beta_ratio"]);
flagged_vessel = awake_beta.vessel_key(awake_beta.beta_ratio < param.flag_beta_range(1) | ...
    awake_beta.beta_ratio > param.flag_beta_range(2));
rows.flagged = ismember(rows.vessel_key, flagged_vessel);
fprintf("flagged %d of %d vessels on the awake fit\n", numel(unique(flagged_vessel)), ...
    numel(unique(rows.vessel_key)));

fig_flagsplit = figure('Color', 'w', 'Position', [230 100 760 520], ...
    'Name', 'fig_fwhm_areagap_beta_flagsplit');
draw_delta(axes(fig_flagsplit), rows, param.beta_column, param.paired_ref, ...
    param.paired_ref_label, param.paired_other, param.paired_other_label, param.fontsize, param.beta_label + ", minus awake", ...
    'split', rows.flagged, 'split_name', ["inside 0..1", "flagged"]);

% which sessions carry the flag, with the path to open them in review_session
flag_rows = rows(rows.flagged & rows.state == "awake", :);
flag_rows = sortrows(flag_rows, "beta_ratio");
for k = 1:height(flag_rows)
    fprintf("   %-14s awake beta %+6.3f | %s\n", string(flag_rows.vessel_key(k)), ...
        flag_rows.beta_ratio(k), string(flag_rows.Directory(k)));
end

%% 5d. Split by which way the vessel moves from awake into NREM
% Two groups, one per sign of the within-vessel NREM difference.
%
% THE NREM BAR IN EACH GROUP IS THE DEFINITION, NOT A RESULT. A group made by the
% sign of that difference is above or below zero in it by construction and its
% p-value means nothing. It is drawn because the group has to be shown for the
% other bars to be read. The question is REM and the transitions: if a vessel that
% rises into NREM also rises there, the direction is a property of the vessel.
[vessel, ~, vessel_of] = unique(rows.vessel_key);
beta_awake = column_mean(rows, vessel_of, numel(vessel), param.beta_column, "awake");
beta_nrem = column_mean(rows, vessel_of, numel(vessel), param.beta_column, "nrem");
has_both = isfinite(beta_awake) & isfinite(beta_nrem);
rises_vessel = vessel(has_both & beta_nrem > beta_awake);
rows.rises = ismember(rows.vessel_key, rises_vessel);
fprintf("into NREM : %d vessels rise, %d fall\n", numel(rises_vessel), ...
    sum(has_both) - numel(rises_vessel));

fig_dirsplit = figure('Color', 'w', 'Position', [260 130 760 520], ...
    'Name', 'fig_fwhm_areagap_beta_dirsplit');
draw_delta(axes(fig_dirsplit), rows(ismember(rows.vessel_key, vessel(has_both)), :), ...
    param.beta_column, param.paired_ref, param.paired_ref_label, ...
    param.paired_other, param.paired_other_label, param.fontsize, ...
    param.beta_label + ", minus awake", ...
    'split', rows.rises(ismember(rows.vessel_key, vessel(has_both))), ...
    'split_name', ["falls into NREM", "rises into NREM"]);

%% 5b. The measured beta against awake, as a difference
% The paired blocks draw both bars and leave the reader to take the difference by
% eye. Here the difference IS the bar, so its height is the effect and its zero is
% the null. Each vessel enters once, as beta(state) - beta(awake) on its own two
% fits, which is the same subtraction tablegeneration_main does for state_summary.
%   a vessel missing either half of a pair is absent from that bar only, so the two n differ
fig_delta = figure('Color', 'w', 'Position', [290 160 520 520], 'Name', 'fig_fwhm_areagap_beta_delta');
draw_delta(axes(fig_delta), rows, param.beta_column, param.paired_ref, ...
    param.paired_ref_label, param.paired_other, param.paired_other_label, ...
    param.fontsize, param.beta_label + ", minus awake");

%% 6. Save, when a block above is worth keeping
% save_figure(fig_scatter, dirs.save_dir, "fig_fwhm_areagap_scatter")
% save_figure(fig_gap,     dirs.save_dir, "fig_fwhm_areagap_um")
% save_figure(fig_frac,    dirs.save_dir, "fig_fwhm_areagap_frac")
% save_figure(fig_state,   dirs.save_dir, "fig_fwhm_areagap_state")
% save_figure(fig_paired,  dirs.save_dir, "fig_fwhm_areagap_paired")
% save_figure(fig_flagsplit, dirs.save_dir, "fig_fwhm_areagap_beta_flagsplit")
% save_figure(fig_dirsplit,  dirs.save_dir, "fig_fwhm_areagap_beta_dirsplit")
% save_figure(fig_delta,   dirs.save_dir, "fig_fwhm_areagap_beta_delta")

%% ---------------------------------------------------------------- helpers
function rgb = state_rgb(clee, state_name)
%STATE_RGB  clee.clist by name, grey for a state it has no colour for.
    if isfield(clee.clist, state_name)
        rgb = clee.clist.(state_name);
    else
        rgb = [0.55 0.55 0.55];
    end
end

function draw_delta(ax, rows, column, ref_state, ref_label, other_state, ...
    other_label, fontsize, y_label, opt)
%DRAW_DELTA  One bar per state, each the within-vessel difference against a reference.
%   The difference IS the bar, so its height is the effect and zero is the null. The
%   paired form -- both bars and a line per vessel -- draws the same number as thirty
%   crossing lines and the eye has to take the difference itself.
% IN   ax          1x1 axes
%      rows        table       the stacked session and state rows
%      column      1x1 str     which column the difference is taken on
%      ref_state   1x1 str     subtracted from every other
%      other_state 1xO str     one bar each
%      fontsize    1x1 double
%      y_label     1x1 str
% opt  split       Nx1 logical two bars per state instead of one, false group first.
%                              [] = one bar per state
%      split_name  1x2 str     what the two groups are called
    arguments
        ax
        rows table
        column      (1,1) string
        ref_state   (1,1) string
        ref_label   (1,1) string
        other_state (1,:) string
        other_label (1,:) string
        fontsize    (1,1) double
        y_label     (1,1) string
        opt.split      (:,1) logical = logical.empty(0,1)
        opt.split_name (1,2) string = ["", ""]
    end
    if isempty(opt.split)
        group_of = {true(height(rows), 1)};
        group_name = "";
    else
        group_of = {~opt.split, opt.split};
        group_name = opt.split_name;
    end
    shade = [0.86 0.90 0.94; 0.72 0.80 0.88];
    offset = ((1:numel(group_of)) - (numel(group_of) + 1) / 2) * 0.30;
    hold(ax, 'on')
    label = strings(1, numel(other_state));
    handle_group = gobjects(1, numel(group_of));
    for g = 1:numel(group_of)
        one = rows(group_of{g}, :);
        [vessel, ~, vessel_of] = unique(one.vessel_key);
        reference = column_mean(one, vessel_of, numel(vessel), column, ref_state);
        for o = 1:numel(other_state)
            other = column_mean(one, vessel_of, numel(vessel), column, other_state(o));
            difference = other - reference;
            difference = difference(isfinite(difference));
            label(o) = other_label(o);
            if isempty(difference)
                continue
            end
            at = o + offset(g);
            mu = mean(difference);
            sd = std(difference);
            ci = 1.96 * sd / sqrt(numel(difference));
            bar_width = 0.5 / numel(group_of);
            handle_group(g) = bar(ax, at, mu, bar_width, 'FaceColor', shade(g, :), ...
                'EdgeColor', [0.3 0.3 0.3]);
            plot(ax, [at at], mu + [-sd sd], '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1);
            plot(ax, [at at], mu + [-ci ci], '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 2.4);
            jitter = (rand(size(difference)) - 0.5) * 0.20 / numel(group_of);
            scatter(ax, at + jitter, difference, 14, [0.2 0.4 0.7], 'filled', ...
                'MarkerFaceAlpha', 0.5);
            fprintf('%-18s %-10s minus %-6s mean %+.4f | CI95 %+.4f .. %+.4f | signrank %.4g | n %d\n', ...
                group_name(min(g, end)), other_state(o), ref_state, mu, mu - ci, ...
                mu + ci, signrank(difference), numel(difference));
        end
    end
    yline(ax, 0, '-', 'Color', [0.4 0.4 0.4]);
    set(ax, 'XTick', 1:numel(other_state), 'XTickLabel', label, 'FontSize', fontsize, ...
        'Box', 'off')
    xlim(ax, [0.4, numel(other_state) + 0.6])
    ylabel(ax, y_label)
    if numel(group_of) > 1
        legend(ax, handle_group, group_name, 'Location', 'best', 'Box', 'off')
    end
    title(ax, "zero = no move from " + ref_label + ", one point per vessel", ...
        'FontWeight', 'normal', 'FontSize', fontsize - 1)
end


function value = column_mean(rows, vessel_of, n_vessel, column, state_name)
%COLUMN_MEAN  One column averaged within vessel, for one state. NaN where absent.
    on_state = rows.state == state_name & isfinite(rows.(column));
    value = accumarray(vessel_of(on_state), rows.(column)(on_state), ...
        [n_vessel 1], @mean, NaN);
end

function draw_paired(ax, rows, column, column_label, ref_state, ref_label, ...
    other_state, other_label, fontsize, y_label, opt)
%DRAW_PAIRED  One reference state against each of the others, within vessel.
%   A vessel missing either half of a pair is dropped from THAT pair only, so the
%   n under each group is its own and the reference bar can differ between groups.
% IN   ax            1x1 axes
%      rows          table       the stacked session and state rows
%      column        1xC str     one group of pairs per entry
%      column_label  1xC str     what to call each on the abscissa
%      ref_state     1x1 str     the left bar of every pair
%      other_state   1xO str     one pair against the reference each
%      fontsize      1x1 double
%      y_label       1x1 str
% opt  null_at       1x1 double  where the reference line goes. 0 for a difference,
%                                1 for a ratio
    arguments
        ax
        rows table
        column       (1,:) string
        column_label (1,:) string
        ref_state    (1,1) string
        ref_label    (1,1) string
        other_state  (1,:) string
        other_label  (1,:) string
        fontsize     (1,1) double
        y_label      (1,1) string
        opt.null_at  (1,1) double = 0
    end
    hold(ax, 'on')
    at = 0;
    tick_at = zeros(1, 0);
    tick_label = strings(1, 0);
    for c = 1:numel(column)
        for o = 1:numel(other_state)
            pair_state = [ref_state, other_state(o)];
            paired_value = pair_by_vessel(rows, column(c), pair_state);
            if isempty(paired_value)
                continue
            end
            for k = 1:2
                value = paired_value(:, k);
                ci = 1.96 * std(value) / sqrt(numel(value));
                bar(ax, at + k, mean(value), 0.5, 'FaceColor', [0.86 0.90 0.94], ...
                    'EdgeColor', [0.3 0.3 0.3], 'BaseValue', opt.null_at);
                plot(ax, [at+k at+k], mean(value) + [-ci ci], '-', ...
                    'Color', [0.1 0.1 0.1], 'LineWidth', 2.4);
            end
            plot(ax, [at+1 at+2], paired_value', '-', 'Color', [0.6 0.6 0.6 0.5], ...
                'LineWidth', 0.7);
            scatter(ax, repmat([at+1 at+2], size(paired_value, 1), 1), paired_value, ...
                14, [0.2 0.4 0.7], 'filled', 'MarkerFaceAlpha', 0.6);
            difference = paired_value(:, 2) - paired_value(:, 1);
            fprintf('%-26s %s %+.4f -> %s %+.4f | within vessel %+.4f | signrank %.4g | n %d\n', ...
                column_label(c), ref_state, mean(paired_value(:, 1)), other_state(o), ...
                mean(paired_value(:, 2)), median(difference), signrank(difference), ...
                size(paired_value, 1));
            tick_at(end+1) = at + 1.5; %#ok<AGROW>
            % the tick names the STATE. With one column its label is the group name,
            % which is a property of the whole figure and belongs in the title -- put
            % on every tick it read "falls into nrem nrem"
            if isscalar(column)
                tick_label(end+1) = other_label(o) + ...
                    " (n=" + size(paired_value, 1) + ")"; %#ok<AGROW>
            else
                tick_label(end+1) = column_label(c) + newline + other_label(o) + ...
                    " (n=" + size(paired_value, 1) + ")"; %#ok<AGROW>
            end
            at = at + 3;
        end
    end
    yline(ax, opt.null_at, '-', 'Color', [0.4 0.4 0.4]);
    % XTickLabel splits a string on newline into separate labels, so the two lines go
    % over as a cellstr with one cell per tick
    set(ax, 'XTick', tick_at, 'XTickLabel', cellstr(split_lines(tick_label)), ...
        'FontSize', fontsize - 1, 'Box', 'off')
    xlim(ax, [0.4, at - 0.4])
    ylabel(ax, y_label)
    heading = "left bar " + ref_label + ", right bar the other, one line per vessel";
    if isscalar(column)
        heading = column_label(1) + "  |  " + heading;
    end
    title(ax, heading, 'FontWeight', 'normal', 'FontSize', fontsize - 1)
end

function stacked = split_lines(label)
%SPLIT_LINES  A two-line tick label as the cellstr MATLAB wants for one tick.
    stacked = cell(1, numel(label));
    for k = 1:numel(label)
        stacked{k} = char(replace(label(k), newline, " "));
    end
    stacked = string(stacked);
end

function paired_value = pair_by_vessel(rows, column, state_pair)
%PAIR_BY_VESSEL  Nx2, one row per vessel that carries both states of the pair.
    on_pair = ismember(rows.state, state_pair);
    one = rows(on_pair, :);
    if isempty(one)
        paired_value = [];
        return
    end
    [vessel, ~, vessel_of] = unique(one.vessel_key);
    paired_value = nan(numel(vessel), 2);
    for k = 1:2
        on_state = one.state == state_pair(k) & isfinite(one.(column));
        if ~any(on_state)
            continue
        end
        paired_value(:, k) = accumarray(vessel_of(on_state), one.(column)(on_state), ...
            [numel(vessel) 1], @mean, NaN);
    end
    paired_value = paired_value(all(isfinite(paired_value), 2), :);
end


function [per_vessel, label] = by_state(rows, column, state_name, state_label)
%BY_STATE  One value per vessel per state and side, sessions averaged within vessel.
%   IN   rows        table       the stacked session and state rows
%        column      1xC str     which columns to pull, one group of bars each
%        state_name  1xS str     the states, in the order they are drawn
%        state_label 1xS str     what each is called on the abscissa
%   OUT  per_vessel  CxS cell    each cell a vector, one entry per vessel
%        label       1xS str     the state and how many vessels it kept
    per_vessel = cell(numel(column), numel(state_name));
    label = strings(1, numel(state_name));
    for s = 1:numel(state_name)
        on_state = rows.state == state_name(s);
        if ~any(on_state)
            label(s) = state_label(s) + " (0)";
            continue
        end
        one_state = rows(on_state, :);
        [~, ~, vessel_of] = unique(one_state.vessel_key);
        for c = 1:numel(column)
            value = accumarray(vessel_of, one_state.(column(c)), [], @mean, NaN);
            per_vessel{c, s} = value(isfinite(value));
        end
        label(s) = state_label(s) + " (n=" + numel(per_vessel{1, s}) + ")";
    end
end

function draw_bars(ax, per_vessel, label, fontsize, y_label, subtitle_text)
%DRAW_BARS  Mean, 95% CI, one SD, and every vessel behind them, C bars per state.
%   The two whiskers say different things and are drawn differently on purpose: the
%   CI is how well the mean is pinned, the SD is how much the vessels differ. A bar
%   with a tight CI over a wide SD is a real mean over a heterogeneous set.
    [n_bar, n_state] = size(per_vessel);
    shade = [0.86 0.90 0.94; 0.72 0.80 0.88];
    offset = ((1:n_bar) - (n_bar + 1) / 2) * 0.30;
    hold(ax, 'on')
    for s = 1:n_state
        for c = 1:n_bar
            value = per_vessel{c, s};
            if isempty(value)
                continue
            end
            at = s + offset(c);
            mu = mean(value);
            sd = std(value);
            ci = 1.96 * sd / sqrt(numel(value));
            bar(ax, at, mu, 0.26, 'FaceColor', shade(min(c, end), :), ...
                'EdgeColor', [0.3 0.3 0.3]);
            plot(ax, [at at], mu + [-sd sd], '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 1);
            plot(ax, [at at], mu + [-ci ci], '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 2.4);
            jitter = (rand(size(value)) - 0.5) * 0.16;
            scatter(ax, at + jitter, value, 14, [0.2 0.4 0.7], 'filled', ...
                'MarkerFaceAlpha', 0.5);
            fprintf('%-14s %-12s mean %+8.3f | CI95 %+8.3f .. %+8.3f | SD %6.3f | n %d\n', ...
                label(s), bar_name(c), mu, mu - ci, mu + ci, sd, numel(value));
        end
    end
    yline(ax, 0, '-', 'Color', [0.4 0.4 0.4]);
    set(ax, 'XTick', 1:n_state, 'XTickLabel', label, 'FontSize', fontsize, ...
        'Box', 'off', 'TickLabelInterpreter', 'none')
    xlim(ax, [0.4, n_state + 0.6])
    ylabel(ax, y_label)
    if subtitle_text ~= ""
        title(ax, subtitle_text, 'FontWeight', 'normal', 'FontSize', fontsize - 1)
    end
end

function name = bar_name(c)
    side = ["constricted", "dilated"];
    name = side(min(c, numel(side)));
end

function save_figure(fig, save_dir, fig_name) %#ok<DEFNU>
% R2023b exportgraphics cannot write svg, so the vector copy goes through print.
    if ~isfolder(save_dir)
        mkdir(save_dir)
    end
    print(fig, fullfile(save_dir, fig_name + ".svg"), '-dsvg', '-vector');
    exportgraphics(fig, fullfile(save_dir, fig_name + ".png"), 'Resolution', 200);
    fprintf('wrote %s\n', fullfile(save_dir, fig_name));
end
