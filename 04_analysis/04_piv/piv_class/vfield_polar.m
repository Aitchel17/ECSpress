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
%     E = (J + J')/2      J = [du/dx du/dy; dv/dx dv/dy]   E is J's symmetric half
%     strain_radial = n' E n        strain_hoop = t' E t        t = n rotated 90
%     strain_shear  = n' E t        rotation = (J - J')/2, one number, no direction
%
%   E rather than J because the antisymmetric half is a rotation, and x'*Om*x = 0
%   for any x -- so a contraction onto a direction cannot see it and building E
%   changes no value here. It is built anyway, because the shear n' E t CAN see it
%   and reads the same three components; taking that one off J would silently
%   subtract the rotation.
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
%   ORDER   applydelaunay -> trifilter -> placetri -> measure -> gate_wedge ->
%           accumulate. Each refuses to run before the one it needs, and they split
%           by what they depend on: applydelaunay on WHERE the vectors are,
%           trifilter on the cells' shapes, placetri on the mask alone, measure on
%           what the vectors SAY. Both filters return a verdict and measure is the
%           one place either is applied.
%           accumulate REQUIRES the verdict as an argument, so a gate cannot be
%           formed and then forgotten.
%
%   see also VFIELD_PROFILE, VFIELD_APPLYDELAUNAY, VFIELD_TRIFILTER,
%            VFIELD_PLACETRI, VFIELD_TRIGRADIENT, SECTOR_WEDGE

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

        mesh  = []     % struct, set by applydelaunay(). See applydelaunay
        keep  = []     % n_tri x 1 logical, set by trifilter(). measure() is the
                       % one place it is applied, so mesh keeps every cell
        place = []     % struct of COLUMNS, set by placetri(). See placetri
        tri   = []     % struct of COLUMNS, set by measure(). See measure
        info = []      % struct; applydelaunay writes the mesh half, measure the
                       % polar half. See both
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

        function applydelaunay(obj)
        %APPLYDELAUNAY  Triangulate the vectors. Every cell, none judged yet.
        %   Sets obj.mesh. Depends on WHERE the vectors are and not on what they
        %   say, so two events over one gate set get the same cells.
        %     mesh.xy   n_point x 2 float px  window centres that carried a vector
        %     mesh.uv   n_point x 2 float px  their displacement, xy's own order
        %     mesh.conn n_tri x 3 int         vertex indices into mesh.xy
        %     mesh.cxy  n_tri x 2 float px    centroid, where the gradient applies
        %     mesh.area n_tri x 1 float px^2  UNSCALED here; measure() converts
        %     mesh.grid_idx n_point x 1 int    linear index into the ny x nx
        %               interrogation grid, so a vertex can be read back against
        %               anything the gates recorded there
        %
        %   note  the cells Delaunay bridges across the masked vessel are IN here.
        %         trifilter takes them out, and they can be seen because this stage
        %         does not
            [obj.mesh, mi] = vfield_applydelaunay(obj.xyuv);
            obj.keep = true(size(obj.mesh.conn, 1), 1);
            obj.info = struct('n_vector', mi.n_vector);
            if obj.param.verbose
                fprintf('%-12s %4d vectors -> %d cells, unjudged\n', ...
                    'delaunay', mi.n_vector, size(obj.mesh.conn, 1));
            end
        end

        function [keep, dropped] = trifilter(obj)
        %TRIFILTER  Which cells are fit to differentiate on. A verdict, not a change.
        %   Sets obj.keep and the dropped counts in obj.info, and returns both so a
        %   caller can look at the verdict without reaching into the object. The
        %   mesh is NOT edited: measure() is the one place the mask is applied, so
        %   a figure can still draw what was thrown away.
        % OUT keep     n_tri x 1 logical
        %     dropped  long_edge / sliver / masked / kept / total, and
        %              max_edge_measured px
            if isempty(obj.mesh)
                error('vfield_polar:notTriangulated', ...
                    'call applydelaunay() first; trifilter() judges its cells.');
            end
            [keep, dropped] = vfield_trifilter(obj.mesh.xy, obj.mesh.conn, ...
                max_edge = obj.param.max_edge, ...
                min_angle = obj.param.min_angle, ...
                exclmask = obj.exclmask);
            obj.keep = keep;
            obj.info.n_tri_kept        = dropped.kept;
            obj.info.max_edge_measured = dropped.max_edge_measured;
            obj.info.dropped           = dropped;
            if obj.param.verbose
                fprintf('%-12s %4d of %d cells kept | dropped %d long, %d sliver, %d masked\n', ...
                    'trifilter', dropped.kept, dropped.total, dropped.long_edge, ...
                    dropped.sliver, dropped.masked);
            end
        end

        function placetri(obj)
        %PLACETRI  Give every triangle its place in the wall's frame. Geometry only.
        %   Sets obj.place. Nothing here reads a displacement, so two events over
        %   one mask get identical output -- which is what lets a claim about a
        %   wedge mean the same thing in both.
        %     r      n_tri x 1 float  MICRONS from the wall
        %     wedge  n_tri x 1 float  1..n_wedge, NaN where the triangle has none
        %     bin    n_tri x 1 float  1..n_bin, NaN outside the bin range
        %     nxy    n_tri x 2 float  [nx ny] outward unit wall normal
            if isempty(obj.mesh)
                error('vfield_polar:notTriangulated', ...
                    'call applydelaunay() first; placetri() places its cells.');
            end
            obj.place = vfield_placetri(obj.mesh.cxy, obj.dist, obj.near_idx, ...
                obj.wedgemap, obj.pixel2um, obj.bin_edges_um);
            if obj.param.verbose
                fprintf('%-12s %4d triangles | %d with a wedge, %d in the bin range\n', ...
                    'placetri', size(obj.mesh.cxy, 1), ...
                    nnz(isfinite(obj.place.wedge)), nnz(isfinite(obj.place.bin)));
            end
        end

        function measure(obj)
        %MEASURE  Differentiate on the mesh and resolve against the wall normal.
        %   The only stage that reads what the vectors SAY. applydelaunay decided
        %   which cells exist and placetri decided where they are; this one decides
        %   what they measured, and it is the only place a physical claim is made.
        %   Sets obj.tri, a struct of COLUMNS all n_tri x 1 (or x 2, x 4) so that
        %   one logical selection applies to every field at once -- which is what
        %   stops a subscript vector and a value vector being indexed differently.
        %     conn         n_tri x 3 int    vertex indices into obj.mesh.xy
        %     divergence   n_tri x 1 float  dimensionless, trace of the gradient
        %     grad         n_tri x 4 float  [du/dx du/dy dv/dx dv/dy], per px
        %     area         n_tri x 1 float  MICRONS^2
        %     cxy          n_tri x 2 float  [x y] px, where it applies
        %     nxy          n_tri x 2 float  copied from obj.place, so one selection
        %                  reaches the geometry and the measurement together
        %     r            n_tri x 1 float  MICRONS from the wall
        %     wedge        n_tri x 1 float  1..n_wedge, NaN outside
        %     bin          n_tri x 1 float  1..n_bin, NaN outside the bin range
        %     disp_radial  n_tri x 1 float  MICRONS, + = outward
        %     disp_tangential n_tri x 1 float MICRONS, + = counter-clockwise.
        %                  display only; no volume claim reads it
        %     strain_radial n_tri x 1 float dimensionless, n' E n
        %     strain_hoop  n_tri x 1 float  dimensionless, t' E t
        %     strain_shear n_tri x 1 float  dimensionless, n' E t
        %     rotation     n_tri x 1 float  dimensionless, (dv/dx - du/dy)/2,
        %                  + = counter-clockwise, the same sense as disp_tangential
        %
        %   note  strain_radial + strain_hoop = divergence identically, for any
        %         unit n -- the cross terms cancel. Quoting both counts one result
        %         twice
        %   note  divergence and rotation are the two columns no choice of n can
        %         change; rotation is the one the other four cannot see at all
        %   note  strain_radial, strain_hoop and strain_shear are the three
        %         components of E in the (n, t) basis, so they carry the whole
        %         tensor. A caller wanting the principal strains reads those
        %         three -- best from accumulate's area-weighted cells, since an
        %         eigenvalue is not linear and must not be averaged itself
            if isempty(obj.place)
                error('vfield_polar:notPlaced', ...
                    'call placetri() first; measure() resolves against its normals.');
            end
            % THE ONE PLACE the filter's verdict is applied. Every column below
            % is subset by it together, so a subscript and a value cannot end up
            % indexed differently
            conn = obj.mesh.conn(obj.keep, :);
            grad = vfield_trigradient(obj.mesh.xy, conn, obj.mesh.uv);   % J, per triangle

            % The displacement AT the centroid. Exact rather than interpolated:
            % the interpolant is linear, so the mean of the three vertices IS the
            % value at the centroid, and it lands where the gradient applies.
            % Six nodal numbers make an affine map: these are its 2 constants and
            % vfield_trigradient returned the other 4
            uvc = (obj.mesh.uv(conn(:,1), :) ...
                 + obj.mesh.uv(conn(:,2), :) ...
                 + obj.mesh.uv(conn(:,3), :)) / 3;

            p2u = obj.pixel2um;
            nx = obj.place.nxy(obj.keep, 1);         % outward wall normal, n
            ny = obj.place.nxy(obj.keep, 2);         % t = (-ny, nx), n turned 90

            % All five projections of ONE tensor, taken together so the identity is
            % visible: the trace, the two stretches where radial + hoop = trace for
            % any unit n, the shear that completes the (n, t) pair, and the rotation
            % that none of the first four can see
            dudx = grad(:,1);                       % J(1,1)
            dudy = grad(:,2);                       % J(1,2)
            dvdx = grad(:,3);                       % J(2,1)
            dvdy = grad(:,4);                       % J(2,2)

            % E = (J + J')/2, the three independent entries. Transposing averages
            % the two off-diagonals into one number and leaves the diagonal alone,
            % which is why only strain_xy does arithmetic
            strain_xx = dudx;                       % E(1,1)
            strain_xy = (dudy + dvdx) / 2;          % E(1,2) = E(2,1). Half of the
                                                    % engineering shear gamma_xy a
                                                    % finite element text prints
            strain_yy = dvdy;                       % E(2,2)

            divergence    = strain_xx + strain_yy;                                            % tr E, first invariant
            strain_radial = nx.*nx.*strain_xx + 2*nx.*ny.*strain_xy + ny.*ny.*strain_yy;      % n' E n
            strain_hoop   = ny.*ny.*strain_xx - 2*nx.*ny.*strain_xy + nx.*nx.*strain_yy;      % t' E t

            % n' E t, the contraction the antisymmetric half reaches. Taking it off
            % the raw gradient instead returns n' E t - rotation, and nothing in the
            % number says a rotation was subtracted -- which is why E is built above
            strain_shear = nx.*ny.*(strain_yy - strain_xx) + (nx.*nx - ny.*ny).*strain_xy;    % n' E t

            % The antisymmetric half, which holds one number. It shares with
            % divergence the property that no n appears, so neither moves if the
            % wall normal is wrong; what makes this one different is that the
            % other four projections cannot see it at all
            rotation = (dvdx - dudy) / 2;           % Om(2,1), = half the curl

            disp_radial = (uvc(:,1).*nx + uvc(:,2).*ny) * p2u;                                % u(centroid) . n
            % the same displacement on the tangential axis. Carried because the
            % display layer had a tangential map before the class existed and
            % dropping a capability silently is worse than keeping a cheap column;
            % nothing in the volume argument reads it
            disp_tangential = (-uvc(:,1).*ny + uvc(:,2).*nx) * p2u;

            obj.tri = struct( ...
                'conn',            conn, ...
                'divergence',      divergence, ...
                'grad',            grad, ...
                'area',            obj.mesh.area(obj.keep) * p2u^2, ...
                'cxy',             obj.mesh.cxy(obj.keep, :), ...
                'nxy',             obj.place.nxy(obj.keep, :), ...
                'r',               obj.place.r(obj.keep), ...
                'wedge',           obj.place.wedge(obj.keep), ...
                'bin',             obj.place.bin(obj.keep), ...
                'disp_radial',     disp_radial, ...
                'disp_tangential', disp_tangential, ...
                'strain_radial',   strain_radial, ...
                'strain_hoop',     strain_hoop, ...
                'strain_shear',    strain_shear, ...
                'rotation',        rotation);

            % info gains the polar half. applydelaunay set the mesh half and it is
            % not rewritten here -- one field, one writer
            placed = isfinite(obj.place.bin(obj.keep)) ...
                     & isfinite(obj.place.wedge(obj.keep)) & isfinite(divergence);
            obj.info.n_tri_placed     = nnz(placed);
            obj.info.wedge_center_rad = obj.wedge_center_measured();

            if obj.param.verbose
                fprintf('%-12s %4d triangles -> %d placed in a bin and a wedge\n', ...
                    'measure', obj.info.n_tri_kept, nnz(placed));
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
        %              strain_shear    dimensionless, n' E t
        %              rotation        dimensionless, + counter-clockwise
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
                'strain_shear', blank, 'rotation', blank, ...
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
                ss = obj.tri.strain_shear(sel);
                rt = obj.tri.rotation(sel);

                a_sum = accumarray(b, a, [nB 1]);
                flux = accumarray(b, a .* dv, [nB 1]);
                dr_sum = accumarray(b, a .* dr, [nB 1]);
                dt_sum = accumarray(b, a .* dt, [nB 1]);
                sr_sum = accumarray(b, a .* sr, [nB 1]);
                sh_sum = accumarray(b, a .* sh, [nB 1]);
                ss_sum = accumarray(b, a .* ss, [nB 1]);
                rt_sum = accumarray(b, a .* rt, [nB 1]);

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
                cells.strain_shear(thick, w) = ss_sum(thick) ./ a_sum(thick);
                cells.rotation(thick, w) = rt_sum(thick) ./ a_sum(thick);

                % The cumulative flux, over the RUN and not the column. Cumulating
                % from bin 1 would put the mask's empty inner bins in front of it,
                % and one leading NaN makes cumsum return NaN for everything after
                % -- the whole wedge silently disappears. Found by the smoke test;
                % see CLAUDE_LOG.md
                cells.volume_out(run, w) = cumsum(flux(run));
            end
        end

        function curve = collapse_wedges(obj, cells)
        %COLLAPSE_WEDGES  The n_bin x n_wedge cells as one radial curve.
        %   accumulate() deliberately stops at n_bin x n_wedge: which wedges to
        %   combine and how is the caller's, not the partition's. This is that
        %   step, done the way the project already treats an unmeasured cell --
        %   NaN, never zero, and nothing invented between measurements.
        %
        %   ORDER IS THE WHOLE THING. Combining each wedge's CUMULATIVE curve is
        %   wrong in both directions: fill a wedge's empty bin with zero and it
        %   claims nothing happened at a radius it never reached, and its frozen
        %   total then rides along in every bin outside; drop the wedge instead
        %   and one interior gap costs its whole curve. Combining the PER-BIN
        %   values first and accumulating after, a gap costs that wedge that bin
        %   and nothing else.
        % IN   cells   the struct accumulate() returned
        % OUT  curve   volume_out    1 x n_bin  MICRONS^2 per micron of vessel,
        %                            cumulative outward. NaN until the first bin
        %                            any wedge measured
        %              divergence    1 x n_bin  dimensionless, area-weighted over
        %                            the wedges present in that bin
        %              n_wedge_bin   1 x n_bin  int, how many wedges that bin had
        %              area          1 x n_bin  MICRONS^2 actually measured there
        %   caution  the per-bin total is scaled by n_wedge / n_wedge_bin, which
        %            says the wedges that were not measured behave like the ones
        %            that were. That is an assumption and it is the only one here;
        %            n_wedge_bin is returned so a reader can see how hard it is
        %            working, and it works hardest at the outer radii -- see
        %            PIV_PLAN.md
            arguments
                obj
                cells struct
            end
            nB = obj.n_bin;

            % 1. What each cell contributed, in microns^2. divergence is an
            %    area-weighted mean, so the volume it stands for is div * area
            contrib = cells.divergence .* cells.area;
            present = ~isnan(contrib);
            n_wedge_bin = sum(present, 2)';

            % 2. Per bin, across the wedges that measured it
            total = sum(contrib, 2, 'omitnan')';
            area  = sum(cells.area, 2, 'omitnan')';
            scaled = total .* (obj.n_wedge ./ max(n_wedge_bin, 1));
            scaled(n_wedge_bin == 0) = NaN;

            div_bin = sum(contrib, 2, 'omitnan')' ./ area;
            div_bin(n_wedge_bin == 0) = NaN;

            % 3. Accumulate outward. A bin nobody measured contributes nothing to
            %    the running total rather than voiding it, but the total stays NaN
            %    until the first bin that WAS measured -- before that there is no
            %    total, as opposed to a total of zero
            running = cumsum(fillmissing(scaled, 'constant', 0));
            running(cumsum(n_wedge_bin > 0) == 0) = NaN;

            curve = struct('volume_out', running, 'divergence', div_bin, ...
                'n_wedge_bin', n_wedge_bin, 'area', area);
            if obj.param.verbose
                fprintf('collapse      %d of %d bins measured | wedges per bin %d..%d\n', ...
                    nnz(n_wedge_bin > 0), nB, min(n_wedge_bin), max(n_wedge_bin));
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

        function map = paint_cells(obj, values)
        %PAINT_CELLS  Spread per-wedge or per-cell values over the frame. DISPLAY ONLY.
        %   Not an engine: nothing downstream reads this, it exists so a number
        %   can be looked at in the place it came from. Hand the result to
        %   showpiv.plot_overlay, which is what draws it.
        %
        %   Not sector_paint. That one takes a FUSED label map and one value per
        %   label; this takes the ring and wedge maps separately and a n_bin x
        %   n_wedge matrix, which is how sector_polar's ordering trap is avoided
        %   -- see its header. The two cannot be merged without choosing the
        %   fused ordering that was deliberately not chosen.
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
