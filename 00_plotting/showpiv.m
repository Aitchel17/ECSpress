classdef showpiv < make_fig
% IN   fig_name  figure name, as make_fig
% OUT  obj       one figure, one axis, PIV layers drawn on it
%
%   make_fig owns the figure, the axis and the house style; this adds the layers
%   a displacement field needs. ONE OBJECT IS ONE PANEL -- a multi-panel figure
%   composes several, rather than the class owning a tiled layout.
%
%   The background goes down as truecolor RGB, so the axis colormap stays free
%   for the scalar overlay. Colours come from color_lee.
%
%   FOUR WAYS TO PUT A SCALAR ON THE IMAGE, and they are not interchangeable
%     plot_overlay  a dense H x W map, one pixel one value. For something that
%                   really is defined per pixel
%     plot_grid     one flat square per interrogation window. The honest display
%                   of PIV output: a vector belongs to a WINDOW, and this is the
%                   only one that draws that window at its true size
%     plot_tri      one flat triangle per Delaunay cell, for vfield_polar
%     plot_contour  interpolates first, so it smooths. Reads well, and it is the
%                   one that can show a boundary that is not in the data
%   err  plot_overlay on a sector-painted map draws the PARTITION as if it were
%        resolution -- flat inside every cell, hard edge at every wall
%
%   plot_pivwall is not one of the four: it draws no value. It is the reference
%   the polar quantities are defined against, and every panel showing one needs it.
%
%   USAGE
%     P = showpiv('ev2 vectors');
%     P.plot_background(S(:,:,706));
%     P.plot_pivwall(coremask);
%     P.plot_quiver(piv.uv);
%     P.plot_overlay(vrmap, 2);          % symmetric limits at +/- 2
%     P.plot_grid(pl.xtable, pl.ytable, div_grid, 0.05);
%     P.plot_quiver(uv_die, 'red', 'red');   % second call, second colour

    properties
        scale   (1,1) double = 5        % arrow length = displacement(px) * scale
        alpha   (1,1) double = 0.55     % scalar overlay transparency
        linew   (1,1) double = 1
        head_len   (1,1) double = 4     % arrowhead length, px on the image
        head_width (1,1) double = 3     % arrowhead base width, px
        head_notch (1,1) double = 0     % 0 = plain triangle, ->1 = swept back barbs
        fill_color (1,:) char = 'magenta'  % head face,  a color_lee.clist field name
        line_color (1,:) char = 'orange'   % shaft and head outline
        cmap                            % overlay colormap, a color_lee gradient
    end

    methods
        function obj = showpiv(fig_name, scale)
        % IN   fig_name  figure name, as make_fig
        %      scale     arrow length multiplier, default 5. Set it here rather
        %                than after the fact: it is what makes two figures
        %                comparable, so it belongs with the figure's identity
            arguments
                fig_name char
                scale (1,1) double {mustBePositive} = 5
            end
            obj@make_fig(fig_name, 'normal');
            obj.scale = scale;
            c = color_lee();
            % Diverging by default: plot_overlay runs -lim..+lim, so the quantity
            % is read against zero and a sequential map would put its bright end
            % on one sign
            obj.cmap = c.gradient.hilo;
        end

        function plot_background(obj, img)
        % IN   img  H x W frame the field is drawn over
        %
        %   Laid down as RGB rather than indexed, so plot_overlay can own the
        %   axis colormap. Contrast is stretched here, not by imcontrast, which
        %   would open a dialog every call.
            g = mat2gray(double(img));
            image(obj.ax, repmat(g, [1 1 3]));
            axis(obj.ax, 'image');
            hold(obj.ax, 'on');
            obj.preset_axis
        end

        function h = plot_pivwall(obj, mask, wall_color)
        % IN   mask        H x W bool  the region whose boundary is drawn
        %      wall_color  char        a color_lee.clist field name
        % OUT  h           the line, one object however many regions there are
        %
        %   The wall is not a measurement, it is the reference the polar
        %   quantities are DEFINED against: radial and hoop are components along
        %   and across it, so a panel carrying either is unreadable without it.
        %   Any mask goes in -- the lumen, the exclusion mask, whatever the panel
        %   is making a claim about.
        %
        %   Traced as a boundary rather than drawn from bwperim. Perimeter pixels
        %   plot as one marker each, which thins to a dashed ring the moment the
        %   axis is not at full resolution; a traced loop stays a closed curve at
        %   any zoom. One line object with NaN between regions, the same way
        %   plot_quiver draws its shafts.
        %
        %   Green by default because obj.cmap runs blue to red, and an outline in
        %   either of those reads as a value.
            arguments
                obj
                mask       (:,:) logical
                wall_color (1,:) char = 'green'
            end
            outlines = bwboundaries(mask);
            if isempty(outlines)
                h = gobjects(0);
                return
            end
            % bwboundaries closes each loop itself, so the NaN separates one
            % region from the next and nothing else
            segments = cell(numel(outlines), 1);
            for k = 1:numel(outlines)
                loop = outlines{k};
                xy_loop = [loop(:,2), loop(:,1)];        % [row col] -> [x y]
                segments{k} = [xy_loop; NaN, NaN];
            end
            xy = vertcat(segments{:});
            edge = color_lee.clist.(wall_color);
            hold(obj.ax, 'on');
            h = line(obj.ax, xy(:,1), xy(:,2), 'Color', edge, 'LineWidth', obj.linew);
            axis(obj.ax, 'image');
        end

        function plot_quiver(obj, uv, fill_color, line_color)
        % IN   uv          H x W x 2 dense map, (:,:,1) = u, (:,:,2) = v, NaN off
        %                  the grid
        %      fill_color  color_lee.clist field name, default obj.fill_color
        %      line_color  the same, for shaft and outline
        %
        %   The two colours are arguments and not just properties so one panel can
        %   carry two populations -- call it once per subset, which is how a gate's
        %   verdict gets drawn (kept in one colour, rejected in another) without
        %   the class knowing anything about gates.
        %
        %   Shaft and head are drawn separately: a NaN-separated line for the
        %   shafts and one patch carrying a filled triangle per arrow. MATLAB's
        %   quiver caps arrows with two open barbs and scales the head with the
        %   arrow, which leaves short vectors headless. Two graphics objects
        %   whatever the vector count.
        %
        %   Nothing is normalised. obj.scale is the only length factor, so the
        %   same setting compares one event against another.
            arguments
                obj
                uv         double
                fill_color (1,:) char = obj.fill_color
                line_color (1,:) char = obj.line_color
            end
            % 0. Vectors that actually landed, at drawing length
            u = uv(:,:,1);  v = uv(:,:,2);
            [r, c] = find(~isnan(u) & ~isnan(v));
            idx = sub2ind(size(u), r, c);
            du = u(idx)*obj.scale;  dv = v(idx)*obj.scale;
            L  = hypot(du, dv);
            in = L > 0;
            r = r(in);  c = c(in);  du = du(in);  dv = dv(in);  L = L(in);

            % 1. Unit direction and its perpendicular
            ex = du./L;   ey = dv./L;
            qx = -ey;     qy =  ex;

            % 2. Head length is absolute and its width fixed, so every arrow reads
            %    the same. A short arrow gives its whole length to the head and
            %    ends up as a bare triangle, which still shows the direction
            hl = min(obj.head_len, L);
            hw = obj.head_width;

            tx = c + du;        ty = r + dv;        % tip
            bx = tx - ex.*hl;   by = ty - ey.*hl;   % head base centre
            nx = bx + ex.*hl*obj.head_notch;        % notch, between base and tip
            ny = by + ey.*hl*obj.head_notch;

            % 3. Shafts, one line object, NaN between arrows. They stop at the notch
            hold(obj.ax, 'on');
            sx = [c, nx, nan(size(c))]';
            sy = [r, ny, nan(size(r))]';
            face = color_lee.clist.(fill_color);
            edge = color_lee.clist.(line_color);
            line(obj.ax, sx(:), sy(:), 'Color', edge, 'LineWidth', obj.linew);

            % 4. Heads, one patch, four vertices per arrow: tip, left barb, notch,
            %    right barb. At notch 0 the middle vertex sits on the base and the
            %    quad renders as a plain triangle. The outline carries the shaft
            %    colour, so the arrow reads as one object over a bright field
            vx = [tx, bx + qx.*hw/2, nx, bx - qx.*hw/2]';
            vy = [ty, by + qy.*hw/2, ny, by - qy.*hw/2]';
            n  = numel(tx);
            patch(obj.ax, 'Faces', reshape(1:4*n, 4, n)', 'Vertices', [vx(:), vy(:)], ...
                  'FaceColor', face, 'EdgeColor', edge, 'LineWidth', obj.linew);
            axis(obj.ax, 'image');
        end

        function [C, h] = plot_contour(obj, x, y, val, levels, smoothing)
        % IN   x, y       scattered sample positions, px
        %      val        the value at each, same length
        %      levels     contour levels, or a count. [] = 8 levels
        %      smoothing  gaussian sigma (px) on the gridded field, 0 = none
        % OUT  C, h       contour matrix and handle, as MATLAB's contour
        %
        %   Interpolates onto the pixel grid first. Contouring a sector-painted
        %   map directly traces the cell edges instead of the field, because that
        %   map is flat inside every cell -- the staircase is the partition, not
        %   the data. Only where samples actually are: the convex hull is filled,
        %   the rest stays NaN so no contour is drawn through empty frame.
            arguments
                obj
                x double
                y double
                val double
                levels = []
                smoothing (1,1) double = 8
            end
            ok = isfinite(x) & isfinite(y) & isfinite(val);
            F  = scatteredInterpolant(x(ok), y(ok), val(ok), 'natural', 'none');
            [H, W] = size(obj.ax.Children(end).CData, 1, 2);
            [X, Y] = meshgrid(1:W, 1:H);
            M = F(X, Y);
            if smoothing > 0
                M = imgaussfilt(M, smoothing, 'FilterDomain', 'spatial');
            end
            hold(obj.ax, 'on');
            if isempty(levels); levels = 8; end
            [C, h] = contour(obj.ax, X, Y, M, levels, 'LineWidth', obj.linew + 0.4);
            colormap(obj.ax, obj.cmap);
            axis(obj.ax, 'image');
        end

        function plot_overlay(obj, M, lim)
        % IN   M    H x W scalar field, NaN where it does not apply
        %      lim  scalar = axis runs -lim..+lim | [lo hi] = as given, for a
        %           quantity with no negative half | [] = symmetric on M
        %
        %   Pass the same lim across events, or each map self-normalises and
        %   amplitude differences between events stop being visible.
            arguments
                obj
                M double
                lim (1,:) double = []
            end
            hold(obj.ax, 'on');
            h = imagesc(obj.ax, M);
            set(h, 'AlphaData', obj.alpha * ~isnan(M));
            colormap(obj.ax, obj.cmap);
            range = showpiv.clim_from(lim, M);
            if ~isempty(range)
                clim(obj.ax, range);
            end
            colorbar(obj.ax);
            axis(obj.ax, 'image');
        end

        function h = plot_grid(obj, x, y, val, lim, cell_px)
        % IN   x, y     window centres in px, any matching shape
        %      val      one value per window, same shape
        %      lim      scalar = -lim..+lim | [lo hi] = as given | [] = symmetric
        %      cell_px  square side in px. [] = the vector spacing, read off x and y
        % OUT  h        the patch, one object however many windows there are
        %
        %   One flat square per interrogation window, drawn at the SIZE OF THE
        %   WINDOW. That is what makes it worth having next to plot_overlay: a PIV
        %   vector is an average over its window, and a display that paints it as
        %   one pixel, or interpolates between windows, claims a resolution the
        %   measurement never had.
        %
        %   caution  the squares are the STEP, not the interrogation area. Windows
        %            overlap 2:1 in this pipeline, so neighbouring squares touch
        %            while the data behind them share half their pixels
            arguments
                obj
                x       double
                y       double
                val     double
                lim     (1,:) double = []
                cell_px = []
            end
            keep = isfinite(x) & isfinite(y) & isfinite(val);
            x = x(keep);
            y = y(keep);
            val = val(keep);
            if isempty(x)
                h = gobjects(0);
                return
            end
            if isempty(cell_px)
                cell_px = showpiv.spacing_of(x, y);
            end

            % 1. Four corners per window, one patch for the lot
            half = cell_px / 2;
            vx = [x(:)-half, x(:)+half, x(:)+half, x(:)-half]';
            vy = [y(:)-half, y(:)-half, y(:)+half, y(:)+half]';
            n  = numel(x);
            hold(obj.ax, 'on');
            h = patch(obj.ax, 'Faces', reshape(1:4*n, 4, n)', ...
                'Vertices', [vx(:), vy(:)], 'FaceVertexCData', val(:), ...
                'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', obj.alpha);
            colormap(obj.ax, obj.cmap);
            range = showpiv.clim_from(lim, val);
            if ~isempty(range)
                clim(obj.ax, range);
            end
            colorbar(obj.ax);
            axis(obj.ax, 'image');
        end

        function h = plot_tri(obj, xy, faces, val, lim)
        % IN   xy     nPoint x 2 [x y] vertex positions, px
        %      faces  nTri x 3 vertex indices -- vfield_polar's obj.tri.conn
        %      val    nTri x 1 value per triangle, or nPoint x 1 per vertex
        %      lim    scalar = -lim..+lim | [lo hi] = as given | [] = symmetric
        % OUT  h      the patch
        %
        %   The display that matches the divergence-theorem route: each triangle
        %   IS the cell whose gradient was computed, so painting it flat draws the
        %   partition the number came from and nothing else. Per-vertex val is
        %   interpolated across the face by MATLAB, which is the only smoothing
        %   anywhere in this class -- pass per-triangle val to avoid it.
            arguments
                obj
                xy    double
                faces double
                val   double
                lim   (1,:) double = []
            end
            % caution  faces first : nTri and nPoint can coincide, and a per-face
            %          value silently interpolated is a smoothing nobody asked for
            if numel(val) == size(faces, 1)
                shading_mode = 'flat';
            elseif numel(val) == size(xy, 1)
                shading_mode = 'interp';
            else
                error('showpiv:triValueSize', ...
                    'val has %d entries; expected %d (per triangle) or %d (per vertex).', ...
                    numel(val), size(faces, 1), size(xy, 1));
            end
            hold(obj.ax, 'on');
            h = patch(obj.ax, 'Faces', faces, 'Vertices', xy, ...
                'FaceVertexCData', val(:), 'FaceColor', shading_mode, ...
                'EdgeColor', 'none', 'FaceAlpha', obj.alpha);
            colormap(obj.ax, obj.cmap);
            range = showpiv.clim_from(lim, val);
            if ~isempty(range)
                clim(obj.ax, range);
            end
            colorbar(obj.ax);
            axis(obj.ax, 'image');
        end
    end

    methods (Static)
        function label = axis_label(name)
        %AXIS_LABEL  What the colour bar or the y axis says for one per-cell quantity.
        %   The single place that knows the six names vfield_polar emits. A name
        %   that is not here THROWS -- a mistyped quantity would otherwise draw an
        %   empty axis, and an empty axis reads as a measurement of zero.
        % IN   name   char, one of the six
        % OUT  label  char, TeX
        %
        %   why  the label carries the SIGN CONVENTION, which is the part a
        %        reader cannot recover from the number. disp_radial is + outward
        %        and disp_tangential is + counter-clockwise, both set in
        %        vfield_polar.measure, and this is where they become visible
        %   note volume_out is cumulative: a bin's value is the integral out to
        %        that radius, not the value at it
            arguments
                name (1,:) char
            end
            % name | axis label. The names are the settled vocabulary; see BACKLOG.md
            known = {'divergence',      'divergence  \DeltaA/A'; ...
                     'strain_radial',   'strain_radial  n''En'; ...
                     'strain_hoop',     'strain_hoop  t''Et'; ...
                     'disp_radial',     'disp_radial  (+ outward, \mum)'; ...
                     'disp_tangential', 'disp_tangential  (+ CCW, \mum)'; ...
                     'volume_out',      'volume_out  (cumulative, \mum^2)'};
            hit = strcmp(known(:, 1), name);
            if ~any(hit)
                error('showpiv:unknownQuantity', '%s is not one of: %s', ...
                    name, strjoin(known(:, 1)', ', '));
            end
            label = known{hit, 2};
        end
    end

    methods (Static, Access = private)
        function range = clim_from(lim, val)
        %CLIM_FROM  Colour axis from the lim argument the three scalar layers take.
        % IN   lim    [] = symmetric on val | scalar = -lim..+lim | [lo hi] = as given
        %      val    the data, read only when lim is empty
        % OUT  range  1 x 2, or [] when there is nothing worth setting
        %
        %   A scalar says the quantity is read against zero, which is why it was
        %   the only form for a long time. A magnitude has no negative half, and
        %   on a symmetric axis it spends half the colormap on values it cannot
        %   take -- so two values are passed straight through.
            if isempty(lim)
                peak = max(abs(val(:)), [], 'omitnan');
                range = [-peak, peak];
            elseif isscalar(lim)
                range = [-lim, lim];
            else
                range = reshape(lim, 1, []);
            end
            if ~all(isfinite(range)) || range(end) <= range(1)
                range = [];
            end
        end

        function step = spacing_of(x, y)
        %SPACING_OF  Vector spacing in px, from the positions alone.
        %   The gap between occupied columns and rows, median of both. Same
        %   quantity vfield_grid_step reads off a sparse map, but from a scattered
        %   list, so the plotting layer does not need the map.
            gap = @(v) median(diff(unique(round(v(:)))));
            step = median([gap(x), gap(y)], 'omitnan');
            if ~isfinite(step) || step < 1; step = 1; end
        end
    end
end
