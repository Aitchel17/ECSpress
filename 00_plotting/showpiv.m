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
%     plot_tri      one flat triangle per Delaunay cell, for vfield_strain('tri')
%     plot_contour  interpolates first, so it smooths. Reads well, and it is the
%                   one that can show a boundary that is not in the data
%   err  plot_overlay on a sector-painted map draws the PARTITION as if it were
%        resolution -- flat inside every cell, hard edge at every wall
%
%   USAGE
%     P = showpiv('ev2 vectors');
%     P.plot_background(S(:,:,706));
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
        %      lim  colour limit; the axis runs -lim..+lim, [] = symmetric on M
        %
        %   Pass the same lim across events, or each map self-normalises and
        %   amplitude differences between events stop being visible.
            arguments
                obj
                M double
                lim = []
            end
            if isempty(lim); lim = max(abs(M(:)), [], 'omitnan'); end
            hold(obj.ax, 'on');
            h = imagesc(obj.ax, M);
            set(h, 'AlphaData', obj.alpha * ~isnan(M));
            colormap(obj.ax, obj.cmap);
            if isfinite(lim) && lim > 0; clim(obj.ax, [-lim lim]); end
            colorbar(obj.ax);
            axis(obj.ax, 'image');
        end

        function h = plot_grid(obj, x, y, val, lim, cell_px)
        % IN   x, y     window centres in px, any matching shape
        %      val      one value per window, same shape
        %      lim      colour limit, axis runs -lim..+lim. [] = symmetric on val
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
                lim     = []
                cell_px = []
            end
            keep = isfinite(x) & isfinite(y) & isfinite(val);
            x = x(keep);  y = y(keep);  val = val(keep);
            if isempty(x); h = gobjects(0); return; end
            if isempty(cell_px); cell_px = obj.spacing_of(x, y); end
            if isempty(lim);     lim = max(abs(val), [], 'omitnan'); end

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
            if isfinite(lim) && lim > 0; clim(obj.ax, [-lim lim]); end
            colorbar(obj.ax);
            axis(obj.ax, 'image');
        end

        function h = plot_tri(obj, xy, faces, val, lim)
        % IN   xy     nPoint x 2 [x y] vertex positions, px
        %      faces  nTri x 3 vertex indices -- vfield_strain('tri') info.tri
        %      val    nTri x 1 value per triangle, or nPoint x 1 per vertex
        %      lim    colour limit, axis runs -lim..+lim. [] = symmetric on val
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
                lim   = []
            end
            if isempty(lim); lim = max(abs(val(:)), [], 'omitnan'); end
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
            if isfinite(lim) && lim > 0; clim(obj.ax, [-lim lim]); end
            colorbar(obj.ax);
            axis(obj.ax, 'image');
        end
    end

    methods (Static, Access = private)
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
