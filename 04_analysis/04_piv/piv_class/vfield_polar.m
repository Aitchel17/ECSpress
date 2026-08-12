classdef vfield_polar < handle
%VFIELD_POLAR  One displacement field, read in polar coordinates about a vessel.
%   Triangulates the vectors once, gives every triangle its radius, its wedge and
%   its radial bin, and resolves its gradient into radial and hoop components
%   against the local WALL NORMAL. accumulate() then bins them. One field, one
%   event; the cohort lives in vfield_profile.
%
%   THE RADIUS IS FROM THE WALL AND THERE IS NO AXIS.
%   bwdist gives distance from the traced boundary, and its second output gives
%   the nearest boundary pixel, so the outward normal n is exact everywhere and
%   costs nothing. Every directional quantity here is a contraction against that
%   n, never a division by a radius:
%
%     strain_radial = n' E n        strain_hoop = t' E t        t = n rotated 90
%
%   The two sum to trace(E) = divergence identically, for ANY unit n -- the cross
%   terms cancel. So there is nothing to check and no residual to store.
%   The axisymmetric form u_r/r would instead need r measured from an AXIS, and a
%   distance-from-the-wall map does not have one; the old code manufactured it as
%   sqrt(nnz(coremask)/pi) and its error is worst at the wall, where the hoop term
%   is largest. See CLAUDE_LOG.md
%
%   RADIAL BINS ARE MICRONS, FIXED BY THE CALLER.
%   Not sector_polar rings. Those take their count from the farthest in-band pixel
%   (an image corner) or their edges from this mask's own sorted distances, so
%   ring k is a different distance in a different recording and rows cannot pool.
%   See CLAUDE_LOG.md
%
%   volume_out IS THE CUMULATIVE FLUX, NOT 2*pi*r*disp_radial.
%   By the divergence theorem the area-weighted sum of divergence over a set of
%   triangles IS the flux through its boundary, so cumulating outward gives the
%   volume that crossed each radius with no geometric assumption at all. The
%   2*pi*r form assumes a circle, and the level set of bwdist around a traced
%   boundary is not one.
%
%   NOTHING IS CUMULATED ACROSS A GAP.
%   A wedge's cumulative runs over the contiguous run of populated bins starting
%   at its first, and is NaN outside it. omitnan here would silently average over
%   a set whose membership changes with radius, which is the error that cost the
%   most on this project. See CLAUDE_LOG.md
%
%   ORDER   measure -> gate_wedge -> accumulate. accumulate REQUIRES the verdict
%           as an argument, so a gate cannot be formed and then forgotten.
%
%   see also VFIELD_PROFILE, VFIELD_GRADIENT_TRI, SECTOR_WEDGE

    properties
        param = struct( ...
            'max_edge',      [], ...   % float px   longest triangle edge kept;
                                ...    %            [] = 3x the median shortest
            'min_angle',     15, ...   % float deg  smallest interior angle kept
            'min_tri_wedge', 10, ...   % int        triangles a wedge needs at all
            'min_tri_bin',   3, ...    % int        triangles a bin needs to count
            'verbose',       true)     % bool       one line per stage. No report
                                       %            method: a line that looks
                                       %            wrong -> obj.info
    end

    properties (SetAccess = private)
        % Settled at construction. Each one decides what the numbers MEAN, so
        % changing it afterwards would not rebuild them.
        % Deliberately not grouped under a container name: the closed vocabulary
        % is param / dirs / info / opt / tmp and these are none of them
        pixel2um       % float             microns per pixel
        n_wedge        % int               angular bins, equal angle
        bin_edges_um   % 1 x nB+1 float    radial bin edges, MICRONS from the wall
        gated          % bool              was xyuv gated before it arrived here

        xyuv           % ny x nx x 4 float (:,:,1:2) = [x y] window centre px,
                       %                   (:,:,3:4) = [u v] displacement px
        coremask       % H x W bool        true inside the vessel
        exclmask       % H x W bool        true = never interpolate here
        dist           % H x W float       px from the core wall
        near_idx       % H x W uint32      linear index of the nearest core pixel
        wedgemap       % H x W float       1..n_wedge, NaN inside the core
        centroid       % 1 x 2 float       [cy cx] px, the mask's centre of mass
        edge_rad       % 1 x n_wedge+1 float  wedge edges, rad

        tri  = []      % struct of COLUMNS, set by measure(). See measure
        info = []      % struct, set by measure(). See measure
    end

    properties (Dependent)
        bin_center_um  % 1 x nB float   midpoint of each radial bin
        n_bin          % int
    end

    methods
        function obj = vfield_polar(xyuv, coremask, pixel2um, opt)
        % IN  xyuv       ny x nx x 4 float px, (:,:,1:2) = [x y] window centre,
        %                (:,:,3:4) = [u v] displacement. ny x nx x N x 4 is summed
        %                over N by the engine. The vectors stay on their own grid;
        %                only the masks below live in frame coordinates
        %     coremask   H x W bool, true inside the vessel
        %     pixel2um   float, microns per pixel
        %     n_wedge    int, angular bins
        %     bin_edges_um  1 x nB+1 float, MICRONS from the wall
        %     exclmask   H x W bool, true = never interpolate
        %     gated      bool, whether xyuv already went through the PIV gates.
        %                Required, and carried into every result: the gates
        %                remove a radius-dependent share of the LARGER vectors,
        %                so a number is not comparable to one taken the other way
            arguments
                xyuv     {mustBeNumeric, mustBeNonempty}
                coremask logical
                pixel2um (1,1) double {mustBePositive}
                opt.n_wedge      (1,1) double {mustBePositive, mustBeInteger} = 12
                opt.bin_edges_um (1,:) double = 0:1.5:40
                opt.exclmask     logical = []
                opt.gated        (1,1) logical
            end
            obj.xyuv = xyuv;
            obj.coremask = coremask;
            obj.pixel2um = pixel2um;
            obj.n_wedge = opt.n_wedge;
            obj.bin_edges_um = opt.bin_edges_um;
            obj.gated = opt.gated;
            if isempty(opt.exclmask)
                obj.exclmask = false(size(coremask));
            else
                obj.exclmask = opt.exclmask;
            end

            [distmap, nearmap] = bwdist(coremask);
            obj.dist = distmap;
            obj.near_idx = nearmap;

            [wmap, winfo] = sector_wedge(coremask, opt.n_wedge);
            obj.wedgemap = wmap;
            obj.centroid = winfo.centroid;
            obj.edge_rad = winfo.wedge_edge_rad;
        end

        function v = get.bin_center_um(obj)
            lo = obj.bin_edges_um(1:end-1);
            hi = obj.bin_edges_um(2:end);
            v = (lo + hi) / 2;
        end

        function v = get.n_bin(obj)
            v = numel(obj.bin_edges_um) - 1;
        end

        function measure(obj)
        %MEASURE  Triangulate, then give every triangle its polar coordinates.
        %   Sets obj.tri, a struct of COLUMNS all n_tri x 1 (or x 2, x 4) so that
        %   one logical selection applies to every field at once -- which is what
        %   stops a subscript vector and a value vector being indexed differently.
        %     divergence   n_tri x 1 float  dimensionless, trace of the gradient
        %     grad         n_tri x 4 float  [du/dx du/dy dv/dx dv/dy], per px
        %     area         n_tri x 1 float  MICRONS^2
        %     cxy          n_tri x 2 float  [x y] px, where it applies
        %     nxy          n_tri x 2 float  [nx ny] outward unit wall normal
        %     r            n_tri x 1 float  MICRONS from the wall
        %     wedge        n_tri x 1 float  1..n_wedge, NaN outside
        %     bin          n_tri x 1 float  1..n_bin, NaN outside the bin range
        %     disp_radial  n_tri x 1 float  MICRONS, + = outward
        %     disp_tangential n_tri x 1 float MICRONS, + = counter-clockwise.
        %                  display only; no volume claim reads it
        %     strain_radial n_tri x 1 float dimensionless, n' E n
        %     strain_hoop  n_tri x 1 float  dimensionless, t' E t
        %   and obj.info: n_vector, n_tri_kept, max_edge_measured, dropped,
        %   wedge_center_rad, n_tri_placed
            [t, gi] = vfield_gradient_tri(obj.xyuv, ...
                max_edge = obj.param.max_edge, ...
                min_angle = obj.param.min_angle, ...
                exclmask = obj.exclmask);

            p2u = obj.pixel2um;
            frame_size = size(obj.dist);

            % where each triangle sits, as a pixel
            px = min(max(round(t.cxy(:,1)), 1), frame_size(2));
            py = min(max(round(t.cxy(:,2)), 1), frame_size(1));
            ci = sub2ind(frame_size, py, px);

            r_px = obj.dist(ci);
            r_um = r_px * p2u;
            wedge = obj.wedgemap(ci);
            bin = discretize(r_um, obj.bin_edges_um);

            % 1. The outward wall normal, from bwdist's own nearest-pixel index.
            %    No axis, no equivalent circle: the direction from the nearest
            %    boundary pixel to here IS the outward normal at the wall
            q = double(obj.near_idx(ci));
            [qy, qx] = ind2sub(frame_size, q);
            dx = t.cxy(:,1) - qx;
            dy = t.cxy(:,2) - qy;
            len = hypot(dx, dy);
            % a centroid landing ON the boundary has no direction. Real, not
            % hypothetical: the vector grid reaches the traced edge
            len(len == 0) = NaN;
            nx = dx ./ len;
            ny = dy ./ len;

            % 2. Resolve the gradient. tangential t is n turned a quarter turn
            dudx = t.grad(:,1);
            dudy = t.grad(:,2);
            dvdx = t.grad(:,3);
            dvdy = t.grad(:,4);
            strain_radial = nx.*nx.*dudx + nx.*ny.*dudy + nx.*ny.*dvdx + ny.*ny.*dvdy;
            strain_hoop   = ny.*ny.*dudx - nx.*ny.*dudy - nx.*ny.*dvdx + nx.*nx.*dvdy;
            disp_radial = (t.uvc(:,1).*nx + t.uvc(:,2).*ny) * p2u;
            % the same displacement on the tangential axis. Carried because the
            % display layer had a tangential map before the class existed and
            % dropping a capability silently is worse than keeping a cheap column;
            % nothing in the volume argument reads it
            disp_tangential = (-t.uvc(:,1).*ny + t.uvc(:,2).*nx) * p2u;

            obj.tri = struct( ...
                'divergence',    t.divergence, ...
                'grad',          t.grad, ...
                'area',          t.area * p2u^2, ...
                'cxy',           t.cxy, ...
                'nxy',           [nx, ny], ...
                'r',             r_um, ...
                'wedge',         wedge, ...
                'bin',           bin, ...
                'disp_radial',     disp_radial, ...
                'disp_tangential', disp_tangential, ...
                'strain_radial',   strain_radial, ...
                'strain_hoop',     strain_hoop);

            placed = isfinite(bin) & isfinite(wedge) & isfinite(t.divergence);
            obj.info = struct( ...
                'n_vector',          gi.n_vector, ...
                'n_tri_kept',        gi.n_tri_kept, ...
                'n_tri_placed',      nnz(placed), ...
                'max_edge_measured', gi.max_edge_measured, ...
                'dropped',           gi.dropped, ...
                'wedge_center_rad',  obj.wedge_center_measured());

            if obj.param.verbose
                fprintf('vfield_polar  %d vectors -> %d tri -> %d placed\n', ...
                    gi.n_vector, gi.n_tri_kept, nnz(placed));
            end
        end

        function [wedge_mask, gate_info] = gate_wedge(obj, opt)
        %GATE_WEDGE  Which wedges carry enough triangles. A verdict, not a change.
        %   Nothing here touches obj.tri. accumulate() takes the verdict as a
        %   required argument, which is the one place it is applied.
        % OUT wedge_mask 1 x n_wedge bool, true = keep
        %     gate_info  n_tri_wedge (1 x n_wedge int), reach_bin (1 x n_wedge),
        %                reason (1 x n_wedge string)
            arguments
                obj
                opt.min_tri_wedge double = []
                opt.bin_range (1,2) double = [NaN NaN]
            end
            if isempty(obj.tri)
                error('vfield_polar:notMeasured', 'call measure() first.');
            end
            floor_n = opt.min_tri_wedge;
            if isempty(floor_n)
                floor_n = obj.param.min_tri_wedge;
            end

            n_tri_wedge = zeros(1, obj.n_wedge);
            first_bin = nan(1, obj.n_wedge);
            reach_bin = nan(1, obj.n_wedge);
            reason = strings(1, obj.n_wedge);
            wedge_mask = false(1, obj.n_wedge);
            for w = 1:obj.n_wedge
                % the same private method accumulate() uses. If these two ever
                % disagreed the verdict would not describe what got binned
                [run_first, run_last, n_here] = obj.wedge_run(w);
                n_tri_wedge(w) = sum(n_here);
                first_bin(w) = run_first;
                reach_bin(w) = run_last;
                if isnan(run_first)
                    reason(w) = "no bin reached the per-bin floor";
                    continue
                end
                if n_tri_wedge(w) < floor_n
                    reason(w) = "too few triangles";
                    continue
                end
                if ~isnan(opt.bin_range(1))
                    if run_first > opt.bin_range(1) || run_last < opt.bin_range(2)
                        reason(w) = "run does not cover the comparison interval";
                        continue
                    end
                end
                reason(w) = "kept";
                wedge_mask(w) = true;
            end

            gate_info = struct('n_tri_wedge', n_tri_wedge, ...
                'first_bin', first_bin, 'reach_bin', reach_bin, 'reason', reason);
            if obj.param.verbose
                fprintf('gate_wedge    %d/%d wedges kept\n', nnz(wedge_mask), obj.n_wedge);
            end
        end

        function cells = accumulate(obj, wedge_mask)
        %ACCUMULATE  Bin the triangles. wedge_mask is REQUIRED, not optional.
        %   Every field is n_bin x n_wedge. A wedge the verdict rejected is all
        %   NaN; so is any bin outside the contiguous populated run, and the
        %   cumulative stops there rather than stepping over it.
        % OUT cells  every field n_bin x n_wedge unless noted:
        %              volume_out      MICRONS^2, cumulative flux outward
        %              divergence      dimensionless, area-weighted mean
        %              disp_radial     MICRONS, + outward
        %              disp_tangential MICRONS, + counter-clockwise. Display only
        %              strain_radial   dimensionless, n' E n
        %              strain_hoop     dimensionless, t' E t
        %              n_tri           int, triangles in the cell
        %              area            MICRONS^2
        %              first_bin       1 x n_wedge, where the run starts
        %              reach_bin       1 x n_wedge, where it ends
            arguments
                obj
                wedge_mask (1,:) logical
            end
            if isempty(obj.tri)
                error('vfield_polar:notMeasured', 'call measure() first.');
            end
            nB = obj.n_bin;
            nW = obj.n_wedge;
            blank = nan(nB, nW);
            cells = struct('volume_out', blank, 'divergence', blank, ...
                'disp_radial', blank, 'disp_tangential', blank, ...
                'strain_radial', blank, 'strain_hoop', blank, ...
                'n_tri', zeros(nB, nW), 'area', blank, ...
                'first_bin', nan(1, nW), 'reach_bin', nan(1, nW));

            for w = 1:nW
                if ~wedge_mask(w)
                    continue
                end
                [run_first, run_last, n_here] = obj.wedge_run(w);
                if isnan(run_first)
                    continue
                end
                % ONE selection, applied to every column. A subscript vector and
                % a value vector indexed differently is how this went wrong before
                sel = obj.tri.wedge == w & isfinite(obj.tri.bin) & isfinite(obj.tri.divergence);
                b = obj.tri.bin(sel);
                a = obj.tri.area(sel);
                dv = obj.tri.divergence(sel);
                dr = obj.tri.disp_radial(sel);
                dt = obj.tri.disp_tangential(sel);
                sr = obj.tri.strain_radial(sel);
                sh = obj.tri.strain_hoop(sel);

                a_sum = accumarray(b, a, [nB 1]);
                flux = accumarray(b, a .* dv, [nB 1]);
                dr_sum = accumarray(b, a .* dr, [nB 1]);
                dt_sum = accumarray(b, a .* dt, [nB 1]);
                sr_sum = accumarray(b, a .* sr, [nB 1]);
                sh_sum = accumarray(b, a .* sh, [nB 1]);

                run = run_first:run_last;
                cells.n_tri(:, w) = n_here;
                cells.first_bin(w) = run_first;
                cells.reach_bin(w) = run_last;

                % the per-bin MEANS, only where the bin holds enough triangles for
                % one to mean anything
                thick = run(n_here(run) >= obj.param.min_tri_bin);
                cells.area(thick, w) = a_sum(thick);
                cells.divergence(thick, w) = flux(thick) ./ a_sum(thick);
                cells.disp_radial(thick, w) = dr_sum(thick) ./ a_sum(thick);
                cells.disp_tangential(thick, w) = dt_sum(thick) ./ a_sum(thick);
                cells.strain_radial(thick, w) = sr_sum(thick) ./ a_sum(thick);
                cells.strain_hoop(thick, w) = sh_sum(thick) ./ a_sum(thick);

                % The cumulative flux, over the RUN and not the column. Cumulating
                % from bin 1 would put the mask's empty inner bins in front of it,
                % and one leading NaN makes cumsum return NaN for everything after
                % -- the whole wedge silently disappears. Found by the smoke test;
                % see CLAUDE_LOG.md
                cells.volume_out(run, w) = cumsum(flux(run));
            end
        end

        function tbl = audit_gate(obj, xyuv_ungated)
        %AUDIT_GATE  What the PIV gates took, by radius, against what they left.
        %   obj.xyuv is the gated field; xyuv_ungated is the same correlation with
        %   the gates switched off. Paired by construction, so the difference is
        %   the gates and nothing else. Both are on the interrogation grid, so the
        %   radius of each window is read out of obj.dist at its centre rather
        %   than by indexing a dense map.
        % OUT tbl  n_bin x 5 table : r_um, n_all, rejected_pct, mag_kept_um,
        %          mag_cut_um. mag_cut > mag_kept means the gates are taking the
        %          LARGER displacements at that radius
            arguments
                obj
                xyuv_ungated {mustBeNumeric}
            end
            kept = ~isnan(obj.xyuv(:,:,3));
            have = ~isnan(xyuv_ungated(:,:,3));
            cut = have & ~kept;
            mag = hypot(xyuv_ungated(:,:,3), xyuv_ungated(:,:,4)) * obj.pixel2um;
            % each window's radius, read at its centre
            frame_size = size(obj.dist);
            wx = min(max(round(xyuv_ungated(:,:,1)), 1), frame_size(2));
            wy = min(max(round(xyuv_ungated(:,:,2)), 1), frame_size(1));
            r_um = obj.dist(sub2ind(frame_size, wy, wx)) * obj.pixel2um;
            bin = discretize(r_um, obj.bin_edges_um);

            nB = obj.n_bin;
            n_all = zeros(nB,1);
            rejected_pct = nan(nB,1);
            mag_kept_um = nan(nB,1);
            mag_cut_um = nan(nB,1);
            for j = 1:nB
                here = bin == j & have;
                n_all(j) = nnz(here);
                if n_all(j) == 0
                    continue
                end
                rejected_pct(j) = 100 * nnz(here & cut) / n_all(j);
                mk = mag(here & ~cut);
                mc = mag(here & cut);
                if ~isempty(mk)
                    mag_kept_um(j) = median(mk);
                end
                if ~isempty(mc)
                    mag_cut_um(j) = median(mc);
                end
            end
            tbl = table(obj.bin_center_um(:), n_all, rejected_pct, mag_kept_um, mag_cut_um, ...
                'VariableNames', {'r_um','n_all','rejected_pct','mag_kept_um','mag_cut_um'});
        end

        function map = paint(obj, values)
        %PAINT  Spread per-wedge or per-cell values over the frame. DISPLAY ONLY.
        %   Not an engine: nothing downstream reads this, it exists so a number
        %   can be looked at in the place it came from.
        % IN  values  1 x n_wedge, or n_bin x n_wedge
        % OUT map     H x W float, NaN where no cell
            arguments
                obj
                values double
            end
            map = nan(size(obj.dist));
            r_um = obj.dist * obj.pixel2um;
            bin = discretize(r_um, obj.bin_edges_um);
            if isvector(values)
                for w = 1:obj.n_wedge
                    here = obj.wedgemap == w & isfinite(bin);
                    map(here) = values(w);
                end
                return
            end
            for w = 1:obj.n_wedge
                for j = 1:obj.n_bin
                    here = obj.wedgemap == w & bin == j;
                    map(here) = values(j, w);
                end
            end
        end
    end

    methods (Access = private)
        function centers = wedge_center_measured(obj)
        % Circular mean of the triangles actually placed in each wedge, about the
        % mask centroid. Not the bin centre: a wedge that runs off the frame keeps
        % only the triangles that fit, and where its measurement sits is not where
        % its bin says it should
            centers = nan(1, obj.n_wedge);
            ang = atan2(obj.tri.cxy(:,2) - obj.centroid(1), ...
                        obj.tri.cxy(:,1) - obj.centroid(2));
            for w = 1:obj.n_wedge
                sel = obj.tri.wedge == w;
                if ~any(sel)
                    continue
                end
                a_here = ang(sel);
                mean_sin = mean(sin(a_here), 'omitnan');
                mean_cos = mean(cos(a_here), 'omitnan');
                centers(w) = mod(atan2(mean_sin, mean_cos), 2*pi);
            end
        end

        function [run_first, run_last, n_here] = wedge_run(obj, w)
        % One wedge's contiguous run of populated bins, and its per-bin counts.
        % Both gate_wedge and accumulate call this, so the verdict and the block
        % cannot describe different bins. A bin under min_tri_bin ends the run:
        % everything past a gap is unreachable, and reporting it as reached is how
        % a set whose membership changes with radius gets averaged
            nB = obj.n_bin;
            sel = obj.tri.wedge == w & isfinite(obj.tri.bin) & isfinite(obj.tri.divergence);
            n_here = zeros(nB, 1);
            run_first = NaN;
            run_last = NaN;
            if ~any(sel)
                return
            end
            n_here = accumarray(obj.tri.bin(sel), 1, [nB 1]);
            % the run is over bins that hold ANY triangle. min_tri_bin gates the
            % per-bin MEANS, not this: a bin with two triangles has an untrustworthy
            % mean but its flux is real, and dropping it out of a CUMULATIVE breaks
            % the integral the divergence theorem is giving us. Only a genuinely
            % empty bin is a gap, because there the flux is unknown
            has = n_here >= 1;
            run_first = find(has, 1);
            if isempty(run_first)
                run_first = NaN;
                return
            end
            run_last = run_first;
            while run_last < nB && has(run_last + 1)
                run_last = run_last + 1;
            end
        end
    end
end
