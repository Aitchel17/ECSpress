function h = vfield_polar_labels(div_cells, info, opt)
%VFIELD_POLAR_LABELS  Overlay per-cell divergence values on a polar map (toggleable).
%   Draws each valid cell's value at its anchor in info.cell_xy and adds a "values"
%   toggle button that shows/hides the whole label group. Call it after the polar
%   map is drawn; the labels go on the given axes.
%
% IN   div_cells  n_rings x n_sectors matrix from vfield_divergence_polar
%      info       its info struct; label anchors come from info.cell_xy
%      ax         target axes, default gca
%      format     sprintf format for the value, default '%.3f'
%      color      label RGB, default black
%      fontsize   label font size, default 8
% OUT  h          struct: text (text handles), button (the toggle uicontrol)

arguments
    div_cells double
    info      struct
    opt.ax          = gca
    opt.format (1,:) char = '%.3f'
    opt.color  double = [0 0 0]
    opt.fontsize (1,1) double = 8
end

% 0. Require the anchors from the polar run
if ~isfield(info, 'cell_xy') || isempty(info.cell_xy)
    error('vfield_polar_labels:noXY', ...
        'info.cell_xy missing -- rerun vfield_divergence_polar (block>1 / valid cells).');
end

% 1. Target axes and cell grid
ax = opt.ax;
xy = info.cell_xy;
[nr, ns] = size(div_cells);
hold(ax, 'on');

% 2. One text label per valid cell
ht = gobjects(0);
for ir = 1:nr
    for is = 1:ns
        % 2.1 Skip cells with no value or no anchor
        if isnan(div_cells(ir, is)) || isnan(xy(ir, is, 1)); continue; end
        % 2.2 Write the value at the cell anchor
        ht(end+1) = text(ax, xy(ir, is, 1), xy(ir, is, 2), ...
            sprintf(opt.format, div_cells(ir, is)), ...
            'Color', opt.color, 'FontSize', opt.fontsize, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'Clipping', 'on', 'Tag', 'polarlabel'); %#ok<AGROW>
    end
end

% 3. Toggle button for the label group
fig = ancestor(ax, 'figure');
btn = uicontrol(fig, 'Style', 'togglebutton', 'String', 'values', 'Value', 1, ...
    'Units', 'normalized', 'Position', [0.005 0.005 0.09 0.05], ...
    'Callback', @(s, ~) local_toggle(s, ht));

% 4. Pack the handles
h = struct('text', ht, 'button', btn);
end

% ---------------------------------------------------------------------------
function local_toggle(s, ht)
%LOCAL_TOGGLE  Show or hide the label group from the button state.
if s.Value; st = 'on'; else; st = 'off'; end
set(ht(isgraphics(ht)), 'Visible', st);
end
