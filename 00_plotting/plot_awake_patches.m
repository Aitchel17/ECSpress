function handle_list = plot_awake_patches(ax, time_table, opt)
%PLOT_AWAKE_PATCHES  The puff x locomotion bout tables as patches on an axes.
%   Caller: analysis_awake_makefig
%
% IN   ax           1x1 axes      drawn into; its hold state is put back, and its y limits too
%                                 unless rows is true, where the strip's caller sets them
%      time_table   1x1 struct    table name -> N x 2 [start end] bouts in seconds, as
%                                 bout_puffmove builds it. Read: puff_go_dur*,
%                                 puff_nogo_dur*, puff_moving_dur*, loco_puff, loco_nopuff,
%                                 loco_restpuff, loco_restnopuff
%      opt.rows     1x1 logical   false = locomotion and puff-touched rests shaded the height of
%                                 the axes, puffs as a bar along its top; true = all seven as
%                                 a strip, y = 2 puff, 1 locomotion and rest, each row 0.8 tall
%      opt.ylim     1x2 float     height of the shading when rows is false; [] = ax.YLim
%      opt.bar_frac 1x1 float     height of the puff bar, as a fraction of that height
% OUT  handle_list  K x 1 patch    one patch per category that has bouts, DisplayName = its
%                                 name, sent behind whatever the axes already holds
    arguments
        ax (1,1) matlab.graphics.axis.Axes
        time_table (1,1) struct
        opt.rows (1,1) logical = false
        opt.ylim (1,:) double = []
        opt.bar_frac (1,1) double {mustBeInRange(opt.bar_frac, 0, 1)} = 0.06
    end

    palette = color_lee.clist;
    category = table( ...
        ["loco_restnopuff"; "loco_restpuff"; "loco_nopuff"; "loco_puff"; "puff_nogo"; "puff_go"; "puff_moving"], ...
        [1; 1; 1; 1; 2; 2; 2], ...
        [false; true; true; true; false; false; false], ...
        'VariableNames', {'name', 'strip_row', 'is_shaded'});

    y_range = opt.ylim;
    if isempty(y_range)
        y_range = ax.YLim;
    end

    was_held = ishold(ax);
    hold(ax, 'on');
    handle_list = gobjects(0, 1);
    for k = 1:height(category)
        name = char(category.name(k));
        y_edge = category_extent(name, category.strip_row(k), category.is_shaded(k), ...
            y_range, opt.rows, opt.bar_frac);
        bouts = gather_bouts(time_table, name);
        if isempty(y_edge) || isempty(bouts)
            continue
        end
        start_s = bouts(:, 1)';
        end_s = bouts(:, 2)';
        x_data = [start_s; end_s; end_s; start_s];
        y_column = [y_edge(1); y_edge(1); y_edge(2); y_edge(2)];
        y_data = repmat(y_column, 1, size(bouts, 1));
        handle_list(end+1, 1) = patch(ax, 'XData', x_data, 'YData', y_data, ...
            'FaceColor', palette.(name), 'EdgeColor', 'none', 'DisplayName', name, ...
            'Tag', ['AwakePatch_' name], 'HitTest', 'off'); %#ok<AGROW>
    end

    if ~was_held
        hold(ax, 'off');
    end
    if ~opt.rows
        ylim(ax, y_range);
    end
    if ~isempty(handle_list)
        % uistack sends them one at a time, so the list goes in reversed to keep the drawn
        % order: puff bars above the shading, everything below the traces
        uistack(flipud(handle_list), 'bottom');
    end
end

function y_edge = category_extent(name, strip_row, is_shaded, y_range, use_rows, bar_frac)
%CATEGORY_EXTENT  [low high] one category fills on the y axis, [] if it is not drawn.
%   Caller: plot_awake_patches
%
% IN   name       1 x n char    category name
%      strip_row  1x1 double    its row in a strip
%      is_shaded  1x1 logical   whether it is shaded behind a trace
%      y_range    1x2 float     height of the axes when it is not a strip
%      use_rows   1x1 logical   true = strip layout
%      bar_frac   1x1 float     puff bar height as a fraction of y_range
% OUT  y_edge     1x2 float | []
    y_edge = [];
    if use_rows
        y_edge = strip_row + [-0.4, 0.4];
    elseif startsWith(name, 'puff_')
        bar_height = bar_frac * (y_range(2) - y_range(1));
        y_edge = [y_range(2) - bar_height, y_range(2)];
    elseif is_shaded
        y_edge = y_range;
    end
end

function bouts = gather_bouts(time_table, name)
%GATHER_BOUTS  Every [start end] row of one category; the per-duration puff tables stacked.
%   Caller: plot_awake_patches
%
% IN   time_table  1x1 struct   table name -> N x 2 bouts
%      name        1 x n char   category name
% OUT  bouts       M x 2 float  start and end, seconds; 0 x 2 if the category is empty
    field_name = string(fieldnames(time_table));
    if startsWith(name, 'puff_')
        is_member = startsWith(field_name, [name '_dur']) & ~endsWith(field_name, '_trans');
    else
        is_member = field_name == name;
    end
    table_list = cellfun(@(f) time_table.(f), cellstr(field_name(is_member)), 'UniformOutput', false);
    bouts = vertcat(zeros(0, 2), table_list{:});
end
