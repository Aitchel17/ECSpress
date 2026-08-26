function ax = plot_sleep_strip(state_data, opt)
%PLOT_SLEEP_STRIP  The sleep score as a band of its own, on a fresh thin axes.
%   The only thing done to the data is the join: sleep_score.mat holds raw scoring
%   bins, not bouts -- every row is binwidth_sec wide and abuts the next, so one
%   continuous NREM stretch arrives as dozens of rows. Drawn as a patch each they
%   seam against one another at strip height and the band reads pale and striped.
%
%   Everything else is plot_sleep_patches', including the colours and the 0.4 that
%   keeps clist.awake from going to a black bar. An empty axes already has
%   YLim [0 1] and that function takes it as given, so nothing here sets it.
%
%   IN   state_data  1x1 struct   <state> -> Nx2 [start end]. A field that is not
%                                 Nx2 numeric is left alone
%   opt  parent      figure | tiledlayout | []   [] opens a strip-shaped figure
%        t_axis      1xT double   pass it when the intervals are INDICES. Two bins
%                                 abut at the same number in seconds and one apart
%                                 in indices, so the gap that still counts differs
%        width       1x1 double   inches, only used when parent is []
%        height      1x1 double   inches, only used when parent is []
%   OUT  ax          1x1 axes     a plain axes with no box, grid, ticks or axis
%                                 lines -- just the band. fusion_figure adopts it
%
%   Example
%     sleepstrip = plot_sleep_strip(sleepscore);
%     fusion_figure([sleepstrip, pax_fig.bv.ax], 'height_each', [0.5 1.8]);
    arguments
        state_data (1,1) struct
        opt.parent = []
        opt.t_axis (1,:) double = []
        opt.width  (1,1) double {mustBePositive} = 7
        opt.height (1,1) double {mustBePositive} = 0.5
    end

    parent = opt.parent;
    if isempty(parent)
        parent = figure('Color', 'w', 'Units', 'inches', 'Name', 'sleep_strip', ...
            'Position', [1 1 opt.width opt.height]);
    end
    ax = axes(parent);
    plot_sleep_patches(ax, join_bins(state_data, opt.t_axis), opt.t_axis);
    % nothing but the band: no box, no grid, no ticks and no axis lines. The
    % abscissa belongs to whatever trace this sits above, which carries it once
    set(ax, 'Box', 'off', 'XTick', [], 'YTick', [], 'XColor', 'none', ...
        'YColor', 'none', 'Layer', 'top');
    grid(ax, 'off')
end

function joined_data = join_bins(state_data, t_axis)
%JOIN_BINS  Each state's rows sorted, then every run of touching rows collapsed.
    gap = 0;
    if ~isempty(t_axis)
        gap = 1;
    end
    joined_data = state_data;
    for f = string(fieldnames(state_data))'
        interval = state_data.(f);
        if ~isnumeric(interval) || isempty(interval) || size(interval, 2) ~= 2
            continue
        end
        interval = sortrows(interval, 1);
        joined = interval;      % nothing touching is the most rows there can be
        n_kept = 1;
        for k = 2:size(interval, 1)
            if interval(k, 1) <= joined(n_kept, 2) + gap
                joined(n_kept, 2) = max(joined(n_kept, 2), interval(k, 2));
            else
                n_kept = n_kept + 1;
                joined(n_kept, :) = interval(k, :);
            end
        end
        joined_data.(f) = joined(1:n_kept, :);
    end
end
