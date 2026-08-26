classdef fusion_figure < handle
%FUSION_FIGURE  A board that finished panels are snapped onto, one at a time.
%   Every makefig in this pipeline opens its own window with one axes in it. This
%   holds an empty grid and copies those axes into named cells, so a composite is
%   built by hand and can be taken apart again -- add the wrong panel, drop it, add
%   another. The source windows are left as they are; this copies out of them.
%
%   copyobj carries the axes, its children, its ticks and its aspect ratio. Three
%   things it does not carry, all handled here:
%     a legend is parented to the FIGURE, not the axes, so the two go over in ONE
%     call or the legend arrives with no entries
%     a colormap set on the source FIGURE is invisible to the copy, so it is read
%     back through colormap(source) and re-applied
%     a colorbar does not come across at all and has to be remade, which is why
%     'colorbar' is an option on add rather than something that just works
%
%   Example
%     board = fusion_figure('grid', [3 2], 'width', 9, 'height', 10);
%     fusion_figure.available
%     board.add(cluster_fig.constricted.ax, 'at', [1 1], 'colorbar', true)
%     board.add(cluster_fig.dilated.ax,     'at', [1 2])
%     board.add(kymo_fig.ax,                'at', [2 1], 'span', [1 2])
%     board.add(pax_fig.eps.ax,             'at', [3 1], 'span', [1 2])
%     board.link_x(1)
%     board.save_figure(dirs.save_dir, "figure3")

    properties (SetAccess = private)
        fig                 % 1x1 figure         the board itself
        layout              % 1x1 TiledChartLayout
        panel               % 1xN struct         at 1x2 | span 1x2 | ax handle |
                            %                    source 1x1 str, where it came from
        param               % 1x1 struct         grid 1x2 | width | height | fontsize
    end

    methods
        function obj = fusion_figure(opt)
        %FUSION_FIGURE  An empty board. Nothing is drawn until add is called.
        %   opt  grid      1x2 double   [rows columns] of cells
        %        width     1x1 double   inches
        %        height    1x1 double   inches
        %        fontsize  1x1 double   one type size over panels drawn at different
        %                               figure sizes. [] leaves each as it came
        %        spacing   1x1 str      'compact' | 'tight' | 'loose' | 'none'
            arguments
                opt.grid     (1,2) double {mustBeInteger, mustBePositive} = [2 2]
                opt.width    (1,1) double {mustBePositive} = 7
                opt.height   (1,1) double {mustBePositive} = 7
                opt.fontsize (1,:) double = []
                opt.spacing  (1,1) string = "compact"
            end
            obj.param = opt;
            obj.panel = struct('at', {}, 'span', {}, 'ax', {}, 'source', {});
            obj.fig = figure('Color', 'w', 'Units', 'inches', 'Name', 'fusion_figure', ...
                'Position', [1 1 opt.width opt.height]);
            obj.layout = tiledlayout(obj.fig, opt.grid(1), opt.grid(2), ...
                'TileSpacing', opt.spacing, 'Padding', opt.spacing);
        end

        function add(obj, source_ax, opt)
        %ADD  One finished axes into one cell.
        %   IN   source_ax  1x1 axes | polaraxes   from any open figure
        %   opt  at         1x2 double   [row column], 1-based
        %        span       1x2 double   how many cells it covers
        %        label      1x1 str      prepended to the panel's own title
        %        colorbar   1x1 logical  remake one, since copyobj drops it
            arguments
                obj
                source_ax   (1,1)
                opt.at      (1,2) double {mustBeInteger, mustBePositive}
                opt.span    (1,2) double {mustBeInteger, mustBePositive} = [1 1]
                opt.label   (1,1) string = ""
                opt.colorbar (1,1) logical = false
            end
            obj.checkcell(opt.at, opt.span);
            copied = obj.adopt(source_ax);
            copied.Layout.Tile = (opt.at(1) - 1) * obj.param.grid(2) + opt.at(2);
            copied.Layout.TileSpan = opt.span;
            if ~isempty(obj.param.fontsize)
                fusion_figure.matchfontsize(copied, obj.param.fontsize);
            end
            if opt.label ~= ""
                own_title = join(string(copied.Title.String), " ");
                title(copied, strtrim(opt.label + "  " + own_title));
            end
            if opt.colorbar
                bar_handle = colorbar(copied);
                source_bar = get(source_ax, 'Colorbar');
                if ~isempty(source_bar)
                    bar_handle.Label.String = source_bar.Label.String;
                end
            end
            obj.panel(end+1) = struct('at', opt.at, 'span', opt.span, 'ax', copied, ...
                'source', string(source_ax.Parent.Name));
            fprintf('added at [%d %d] span [%d %d], %d panels on the board\n', ...
                opt.at, opt.span, numel(obj.panel));
        end

        function drop(obj, at)
        %DROP  Take the panel in this cell back off. The source window is untouched.
        %   IN   at  1x2 double   [row column] as it was given to add
            arguments
                obj
                at (1,2) double {mustBeInteger, mustBePositive}
            end
            occupied = arrayfun(@(p) isequal(p.at, at), obj.panel);
            if ~any(occupied)
                fprintf('nothing at [%d %d]\n', at);
                return
            end
            delete(obj.panel(occupied).ax);
            obj.panel(occupied) = [];
            fprintf('dropped [%d %d], %d panels left\n', at, numel(obj.panel));
        end

        function link_x(obj, column)
        %LINK_X  One abscissa down one column, spelled out on the bottom panel only.
        %   Per column and not over the whole board: a board holds panels that share
        %   time beside panels that have no abscissa at all, and linking a polaraxes
        %   to a trace is not a thing that can be done.
        %   IN   column  1x1 double   which column of the grid
            arguments
                obj
                column (1,1) double {mustBeInteger, mustBePositive}
            end
            in_column = arrayfun(@(p) p.at(2) == column && ...
                isa(p.ax, 'matlab.graphics.axis.Axes'), obj.panel);
            if sum(in_column) < 2
                fprintf('column %d has %d linkable panels, nothing to link\n', ...
                    column, sum(in_column));
                return
            end
            stacked = obj.panel(in_column);
            [~, top_first] = sort(arrayfun(@(p) p.at(1), stacked));
            stacked = stacked(top_first);
            linkaxes([stacked.ax], 'x');
            for k = 1:numel(stacked) - 1
                set(stacked(k).ax, 'XTickLabel', []);
                xlabel(stacked(k).ax, '');
            end
            fprintf('column %d : %d panels on one abscissa\n', column, numel(stacked));
        end

        function save_figure(obj, save_dir, fig_name)
        %SAVE_FIGURE  The board as svg and png.
        %   R2023b exportgraphics cannot write svg, so the vector copy goes through
        %   print, which is what make_fig.save2svg already does. see CLAUDE.md
            arguments
                obj
                save_dir (1,1) string
                fig_name (1,1) string
            end
            if ~isfolder(save_dir)
                mkdir(save_dir)
            end
            print(obj.fig, fullfile(save_dir, fig_name + ".svg"), '-dsvg', '-vector');
            exportgraphics(obj.fig, fullfile(save_dir, fig_name + ".png"), 'Resolution', 200);
            fprintf('wrote %s.svg and .png to %s\n', fig_name, save_dir);
        end
    end

    methods (Access = private)
        function checkcell(obj, at, span)
        % the grid is fixed at construction, so a cell outside it is a typo the
        % caller wants told about now and not as an empty tile later
            last = at + span - 1;
            if any(last > obj.param.grid)
                % error takes scalars where fprintf would flatten an array
                error('fusion_figure:offGrid', ...
                    '[%d %d] span [%d %d] reaches cell [%d %d], past the %dx%d grid', ...
                    at(1), at(2), span(1), span(2), last(1), last(2), ...
                    obj.param.grid(1), obj.param.grid(2));
            end
            for k = 1:numel(obj.panel)
                taken = obj.panel(k).at;
                if isequal(taken, at)
                    error('fusion_figure:cellTaken', ...
                        '[%d %d] already holds a panel from %s. drop it first', ...
                        at(1), at(2), obj.panel(k).source);
                end
            end
        end

        function copied = adopt(obj, source)
        % the legend has to travel in the same copyobj call or it arrives empty, and
        % colormap(ax) answers with the FIGURE's map when the axes has none of its
        % own, which is the case a spectrogram falls into
            source_legend = get(source, 'Legend');
            if isempty(source_legend)
                handles = copyobj(source, obj.layout);
            else
                handles = copyobj([source, source_legend], obj.layout);
            end
            copied = handles(1);
            colormap(copied, colormap(source));
            if isprop(source, 'CLim')
                copied.CLim = source.CLim;
            end
        end
    end

    methods (Static)
        function available()
        %AVAILABLE  Every axes on every open figure, with what add would need.
        %   The board cannot list what it can hold, only what is open. Run this
        %   after the makefig calls and before building.
            open_figures = findall(groot, 'Type', 'figure');
            if isempty(open_figures)
                fprintf('no figures open\n');
                return
            end
            fprintf('%-24s %-10s %s\n', 'figure', 'kind', 'title');
            for k = numel(open_figures):-1:1
                one = open_figures(k);
                if isa(one, 'matlab.ui.Figure') && strcmp(one.Name, 'fusion_figure')
                    continue
                end
                found = findall(one, 'Type', 'axes', '-or', 'Type', 'polaraxes');
                for a = 1:numel(found)
                    fprintf('%-24s %-10s %s\n', fusion_figure.shortname(one), ...
                        fusion_figure.kindof(found(a)), ...
                        join(string(found(a).Title.String), " "));
                end
            end
        end
    end

    methods (Static, Access = private)
        function matchfontsize(panel, base_size)
        %MATCHFONTSIZE  One tick size, with the label and title kept in the proportion that
        %   axes was already using rather than flattened to the ticks. A legend carries its
        %   own size and does not follow the axes, so it is set too.
            panel.FontSize = base_size;
            if isprop(panel, 'LabelFontSizeMultiplier')
                panel.XLabel.FontSize = base_size * panel.LabelFontSizeMultiplier;
                panel.YLabel.FontSize = base_size * panel.LabelFontSizeMultiplier;
            end
            panel.Title.FontSize = base_size * panel.TitleFontSizeMultiplier;
            panel_legend = get(panel, 'Legend');
            if ~isempty(panel_legend)
                panel_legend.FontSize = base_size;
            end
        end

        function name = shortname(one)
            name = string(one.Name);
            if name == ""
                name = "Figure " + string(one.Number);
            end
        end

        function kind = kindof(one)
            if isa(one, 'matlab.graphics.axis.PolarAxes')
                kind = "polar";
            elseif ~isempty(findall(one, 'Type', 'image'))
                kind = "image";
            else
                kind = "axes";
            end
        end
    end
end
