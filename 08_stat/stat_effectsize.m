%% STAT_EFFECTSIZE  Effect sizes for state comparisons from statsummary_table.
%   For each metric (DataType) and each pair of sleep states, computes
%   Cohen's d, Hedges' g (bias-corrected) and Cliff's delta (non-parametric),
%   each with a 95% CI, and writes an `effectsize` table (.mat + .csv).
%
%   Uses the `statsummary_table` currently in the workspace (or loads it).
%   Requires Statistics & Machine Learning Toolbox (meanEffectSize, R2022a+).

%% ---------------- config ----------------
norm_field = 'abs_numeric';   % 'abs_numeric' | 'awakenorm_numeric' | 'awakesubtract_numeric'
value_col  = 'raw_mean';      % 'raw_mean' | 'raw_median' | 'Q2Q3_mean'
agg_level  = 'bout';          % 'bout' (all obs) | 'vessel' | 'mouse'  (biological unit)
ref_state  = '';              % '' = all pairwise; else each state vs this one (e.g. 'awake')
save_dir   = fileparts(mfilename('fullpath'));   % 08_stat

%% ---------------- load ----------------
if ~exist('statsummary_table', 'var')
    error('stat_effectsize:noData', 'statsummary_table not found in the workspace.');
end
ft = statsummary_table.(norm_field).filtered_table;

%% ---------------- aggregate to the chosen unit ----------------
switch agg_level
    case 'bout'
        T = ft(:, {'DataType', 'state_name', value_col});
        T = renamevars(T, value_col, 'val');
    case 'vessel'
        T = groupsummary(ft, {'DataType', 'state_name', 'VesselID'}, 'mean', value_col);
        T = renamevars(T, ['mean_' value_col], 'val');
    case 'mouse'
        T = groupsummary(ft, {'DataType', 'state_name', 'MouseID'}, 'mean', value_col);
        T = renamevars(T, ['mean_' value_col], 'val');
    otherwise
        error('stat_effectsize:aggLevel', 'agg_level must be bout | vessel | mouse.');
end

%% ---------------- comparison pairs ----------------
metrics = unique(T.DataType);
states  = unique(T.state_name);
if isempty(ref_state)
    pairs = nchoosek(1:numel(states), 2);                 % all pairwise (A vs B)
else
    ri     = find(states == ref_state, 1);
    others = setdiff(1:numel(states), ri);
    pairs  = [others(:), repmat(ri, numel(others), 1)];   % each state vs the reference (B = ref)
end

%% ---------------- effect sizes ----------------
rows = cell(0, 14);
for mi = 1:numel(metrics)
    dm = metrics(mi);
    for pj = 1:size(pairs, 1)
        sA = states(pairs(pj, 1));
        sB = states(pairs(pj, 2));
        x = T.val(T.DataType == dm & T.state_name == sA);
        y = T.val(T.DataType == dm & T.state_name == sB);
        x = x(~isnan(x));  y = y(~isnan(y));
        nA = numel(x);     nB = numel(y);
        if nA < 2 || nB < 2, continue; end

        ec = meanEffectSize(x, y, Effect = "cohen");   % Cohen's d (+ exact CI)
        d  = ec.Effect;  dlo = ec.ConfidenceIntervals(1);  dhi = ec.ConfidenceIntervals(2);

        J  = 1 - 3 / (4 * (nA + nB) - 9);              % small-sample bias correction
        g  = d * J;                                    % Hedges' g

        ef = meanEffectSize(x, y, Effect = "cliff");   % Cliff's delta (+ bootstrap CI)
        cd = ef.Effect;  clo = ef.ConfidenceIntervals(1);  chi = ef.ConfidenceIntervals(2);

        rows(end + 1, :) = {char(dm), char(sA), char(sB), nA, nB, mean(x), mean(y), ...
            d, dlo, dhi, g, cd, clo, chi}; %#ok<SAGROW>
    end
end

effectsize = cell2table(rows, 'VariableNames', ...
    {'DataType', 'stateA', 'stateB', 'nA', 'nB', 'meanA', 'meanB', ...
     'cohens_d', 'd_CI_low', 'd_CI_high', 'hedges_g', 'cliff_delta', 'cliff_CI_low', 'cliff_CI_high'});
effectsize = sortrows(effectsize, {'DataType', 'cohens_d'}, {'ascend', 'descend'});

%% ---------------- report + save ----------------
disp(effectsize);
tag     = sprintf('%s_%s_%s', norm_field, value_col, agg_level);
outfile = fullfile(save_dir, ['effectsize_' tag '.mat']);
save(outfile, 'effectsize');
writetable(effectsize, replace(outfile, '.mat', '.csv'));
fprintf('[stat_effectsize] %d comparisons -> %s (+ .csv)\n', height(effectsize), outfile);
