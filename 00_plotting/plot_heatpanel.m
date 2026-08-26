function handle = plot_heatpanel(ax, x_value, y_value, counts, opt)
%PLOT_HEATPANEL  One 2-D heatmap on one axes: the map, a curve over it, and
%   reference curves through the origin.
%   Draws and returns handles. It writes no label, no title and no legend, because
%   the three panels that call it say different things about the same picture and
%   the caller is the one that knows which.
%
%   IN   ax          1x1 axes
%        x_value     1xN double    bin centres, already in the unit to be shown
%        y_value     1xM double
%        counts      MxN double    NaN where the column was not filled
%   opt  curve       1xN double    the column statistic, drawn over the map. [] = none
%        ref_curve   RxN double    one reference per ROW, on the same x_value.
%                                  Curves and not slopes: the area-conserving
%                                  reference is a level curve, and its slope runs
%                                  0.42 to 0.65 across the diameters one panel
%                                  spans. A straight line is that curve's tangent
%                                  at the origin stretched over the whole range
%        ref_style   1xR str       one line style per reference
%        ref_color   Rx3 double    one colour per reference
%        colormap    Kx3 double
%        curve_color 1x3 double
%        curve_width 1x1 double
%        fontsize    1x1 double
%        ticks       1x1 logical   false strips both axes, for a tile in a grid
%   OUT  handle      struct        .surface .curve .ref, so the caller can legend them
    arguments
        ax
        x_value     (1,:) double
        y_value     (1,:) double
        counts      (:,:) double
        opt.curve       (1,:) double = []
        opt.ref_curve   (:,:) double = []
        opt.ref_style   (1,:) string = "--"
        opt.ref_color   (:,3) double = []
        opt.colormap    (:,3) double = parula(256)
        opt.curve_color (1,3) double = [1 1 1]
        opt.curve_width (1,1) double = 2
        opt.fontsize    (1,1) double = 11
        opt.ticks       (1,1) logical = true
    end

    % log of the count, with a floor rather than a guard: an empty bin is NaN by
    % the time it arrives here and log10 of it stays NaN, which pcolor leaves as
    % the axes background
    pcolor(ax, x_value, y_value, log10(counts + 1e-4));
    shading(ax, 'flat');
    colormap(ax, opt.colormap);
    hold(ax, 'on')

    handle.curve = gobjects(0);
    if ~isempty(opt.curve)
        handle.curve = plot(ax, x_value, opt.curve, '-', 'Color', opt.curve_color, ...
            'LineWidth', opt.curve_width);
    end

    n_ref = size(opt.ref_curve, 1);
    ref_style = opt.ref_style;
    if isscalar(ref_style)
        ref_style = repmat(ref_style, 1, n_ref);
    end
    ref_color = opt.ref_color;
    if isempty(ref_color)
        ref_color = [1 1 1] .* (0.30 + 0.16 * (1:n_ref))';
    end
    handle.ref = gobjects(1, n_ref);
    for r = 1:n_ref
        handle.ref(r) = plot(ax, x_value, opt.ref_curve(r, :), ref_style(r), ...
            'LineWidth', 1.2, 'Color', ref_color(r, :));
    end

    % black ground, so the colormap's dark end reads as empty rather than as low
    set(ax, 'Color', [0 0 0], 'FontSize', opt.fontsize, 'Layer', 'top', 'Box', 'on')
    if ~opt.ticks
        set(ax, 'XTick', [], 'YTick', [])
    end
    ax.Toolbar.Visible = 'off';   % it exports into the panel otherwise
    handle.surface = ax.Children(end);
end
