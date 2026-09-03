classdef vfield_polar < handle
%VFIELD_POLAR  One displacement field, read in polar coordinates about a vessel.
%   Triangulates the vectors once, gives every triangle its radius (bwdist from the
%   traced wall; there is no axis), its wedge and its radial bin, and resolves its
%   gradient against the local outward wall normal n. accumulate() bins the result.
%   One field, one event; the cohort lives in vfield_profile.
%
%     E = (J + J')/2      J = [du/dx du/dy; dv/dx dv/dy]      t = n rotated 90
%     strain_radial = n' E n    strain_hoop = t' E t    strain_shear = n' E t
%     rotation = (J - J')/2     strain_radial + strain_hoop = trace(E) = divergence
%
%   Radial bins are MICRONS from the wall, fixed by the caller. volume_out is the
%   cumulative flux by the divergence theorem, never 2*pi*r*disp_radial, and nothing
%   is cumulated across a gap. The reasons: see CLAUDE_LOG.md, vfield_polar
%
%   ORDER   applydelaunay -> trifilter -> placetri -> measure -> gate_wedge -> accumulate.
%           Both filters return a verdict: measure() applies trifilter's, and
%           accumulate() takes gate_wedge's as a required argument
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
            'verbose',       true)     % bool       one line per stage
    end

    properties (SetAccess = private)
        % settled at construction; each one decides what the numbers mean
        pixel2um       % float             microns per pixel
        n_wedge        % int               angular bins, equal angle
        bin_edges_um   % 1 x nB+1 float    radial bin edges, MICRONS from the wall
        gate_name      % string            the gates xyuv came through before it arrived here

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
        %                (:,:,3:4) = [u v] displacement, on the interrogation grid
        %     coremask   H x W bool, true inside the vessel
        %     pixel2um   float, microns per pixel
        %     n_wedge    int, angular bins
        %     bin_edges_um  1 x nB+1 float, MICRONS from the wall
        %     exclmask   H x W bool, true = never interpolate
        %     gate_name  string, the PIV gates xyuv came through, joined by +.
        %                Required; carried into every result
            arguments
                xyuv     {mustBeNumeric, mustBeNonempty}
                coremask logical
                pixel2um (1,1) double {mustBePositive}
                opt.n_wedge      (1,1) double {mustBePositive, mustBeInteger} = 12
                opt.bin_edges_um (1,:) double = 0:1.5:40
                opt.exclmask     logical = []
                opt.gate_name    (1,1) string
            end
            obj.xyuv = xyuv;
            obj.coremask = coremask;
            obj.pixel2um = pixel2um;
            obj.n_wedge = opt.n_wedge;
            obj.bin_edges_um = opt.bin_edges_um;
            obj.gate_name = opt.gate_name;
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
        %   Sets obj.mesh from where the vectors are, not from what they say.
        %     mesh.xy   n_point x 2 float px  window centres that carried a vector
        %     mesh.uv   n_point x 2 float px  their displacement, xy's own order
        %     mesh.conn n_tri x 3 int         vertex indices into mesh.xy
        %     mesh.cxy  n_tri x 2 float px    centroid, where the gradient applies
        %     mesh.area n_tri x 1 float px^2  UNSCALED here; measure() converts
        %     mesh.grid_idx n_point x 1 int    linear index into the ny x nx
        %               interrogation grid, so a vertex can be read back against
        %               anything the gates recorded there
        %
        %   note  cells bridging the masked vessel are still in; trifilter takes them out
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
        %   Sets obj.keep and obj.info.dropped; the mesh is not edited, measure() applies it.
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
        %   Sets obj.place; nothing here reads a displacement.
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
        %   The only stage that reads what the vectors say. Sets obj.tri, a struct of
        %   COLUMNS, all n_tri x 1 (or x 2, x 4), so one selection applies to every field.
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
        %     disp_tangential n_tri x 1 float MICRONS, + = counter-clockwise. Display only
        %     strain_radial n_tri x 1 float dimensionless, n' E n
        %     strain_hoop  n_tri x 1 float  dimensionless, t' E t
        %     strain_shear n_tri x 1 float  dimensionless, n' E t
        %     rotation     n_tri x 1 float  dimensionless, (dv/dx - du/dy)/2,
        %                  + = counter-clockwise, the same sense as disp_tangential
        %
        %   note  strain_radial, strain_hoop and strain_shear are E in the (n, t) basis,
        %         the whole tensor. Principal strains: from accumulate's cells, not per
        %         triangle, since an eigenvalue must not be averaged itself
            if isempty(obj.place)
                error('vfield_polar:notPlaced', ...
                    'call placetri() first; measure() resolves against its normals.');
            end
            % the one place trifilter's verdict is applied; every column is subset by it together
            conn = obj.mesh.conn(obj.keep, :);
            grad = vfield_trigradient(obj.mesh.xy, conn, obj.mesh.uv);   % J, per triangle

            % displacement at the centroid: the interpolant is linear, so the vertex mean is exact
            uvc = (obj.mesh.uv(conn(:,1), :) ...
                 + obj.mesh.uv(conn(:,2), :) ...
                 + obj.mesh.uv(conn(:,3), :)) / 3;

            p2u = obj.pixel2um;
            nx = obj.place.nxy(obj.keep, 1);         % outward wall normal, n
            ny = obj.place.nxy(obj.keep, 2);         % t = (-ny, nx), n turned 90

            % the five projections of one tensor; the identities are in the header
            dudx = grad(:,1);                       % J(1,1)
            dudy = grad(:,2);                       % J(1,2)
            dvdx = grad(:,3);                       % J(2,1)
            dvdy = grad(:,4);                       % J(2,2)

            % E = (J + J')/2, the three independent entries
            strain_xx = dudx;                       % E(1,1)
            strain_xy = (dudy + dvdx) / 2;          % E(1,2) = E(2,1), half the engineering shear
            strain_yy = dvdy;                       % E(2,2)

            divergence    = strain_xx + strain_yy;                                            % tr E, first invariant
            strain_radial = nx.*nx.*strain_xx + 2*nx.*ny.*strain_xy + ny.*ny.*strain_yy;      % n' E n
            strain_hoop   = ny.*ny.*strain_xx - 2*nx.*ny.*strain_xy + nx.*nx.*strain_yy;      % t' E t

            % n' E t off E, not off J: off J it would carry -rotation
            strain_shear = nx.*ny.*(strain_yy - strain_xx) + (nx.*nx - ny.*ny).*strain_xy;    % n' E t

            % the antisymmetric half, one number; the four projections above cannot see it
            rotation = (dvdx - dudy) / 2;           % Om(2,1), = half the curl

            disp_radial = (uvc(:,1).*nx + uvc(:,2).*ny) * p2u;                                % u(centroid) . n
            % the same displacement on the tangential axis; display only
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

            % the polar half of info; applydelaunay wrote the mesh half
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
                % the same run accumulate() bins, so verdict and block agree
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
                % one selection, applied to every column
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

                % the per-bin means, only where the bin holds enough triangles
                thick = run(n_here(run) >= obj.param.min_tri_bin);
                cells.area(thick, w) = a_sum(thick);
                cells.divergence(thick, w) = flux(thick) ./ a_sum(thick);
                cells.disp_radial(thick, w) = dr_sum(thick) ./ a_sum(thick);
                cells.disp_tangential(thick, w) = dt_sum(thick) ./ a_sum(thick);
                cells.strain_radial(thick, w) = sr_sum(thick) ./ a_sum(thick);
                cells.strain_hoop(thick, w) = sh_sum(thick) ./ a_sum(thick);
                cells.strain_shear(thick, w) = ss_sum(thick) ./ a_sum(thick);
                cells.rotation(thick, w) = rt_sum(thick) ./ a_sum(thick);

                % the cumulative flux over the run, not the column: a leading NaN would
                % void the whole wedge; see CLAUDE_LOG.md
                cells.volume_out(run, w) = cumsum(flux(run));
            end
        end

        function curve = collapse_wedges(obj, cells)
        %COLLAPSE_WEDGES  The n_bin x n_wedge cells as one radial curve.
        %   Combines the wedges per bin FIRST and accumulates after, so a gap costs
        %   that wedge that bin and nothing else. An unmeasured cell is NaN, never
        %   zero. Why not the other order: see CLAUDE_LOG.md, vfield_polar
        % IN   cells   the struct accumulate() returned
        % OUT  curve   volume_out    1 x n_bin  MICRONS^2 per micron of vessel,
        %                            cumulative outward. NaN until the first bin
        %                            any wedge measured
        %              divergence    1 x n_bin  dimensionless, area-weighted over
        %                            the wedges present in that bin
        %              n_wedge_bin   1 x n_bin  int, how many wedges that bin had
        %              area          1 x n_bin  MICRONS^2 actually measured there
        %   caution  the per-bin total is scaled by n_wedge / n_wedge_bin : unmeasured
        %            wedges are assumed to behave like measured ones. n_wedge_bin is
        %            returned so that assumption is visible
            arguments
                obj
                cells struct
            end
            nB = obj.n_bin;

            % 1. what each cell contributed, microns^2 : divergence is an area-weighted mean
            contrib = cells.divergence .* cells.area;
            present = ~isnan(contrib);
            n_wedge_bin = sum(present, 2)';

            % 2. per bin, across the wedges that measured it
            total = sum(contrib, 2, 'omitnan')';
            area  = sum(cells.area, 2, 'omitnan')';
            scaled = total .* (obj.n_wedge ./ max(n_wedge_bin, 1));
            scaled(n_wedge_bin == 0) = NaN;

            div_bin = sum(contrib, 2, 'omitnan')' ./ area;
            div_bin(n_wedge_bin == 0) = NaN;

            % 3. accumulate outward; an unmeasured bin adds nothing, and the total is
            %    NaN until the first measured bin
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
        %   obj.xyuv is the gated field, xyuv_ungated the same correlation before the
        %   gates; each window's radius is read off obj.dist at its centre.
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
        %   showpiv.plot_overlay draws the result. Not sector_paint: that takes a fused
        %   label map, this takes the ring and wedge maps separately.
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
        % circular mean of the triangles actually placed in each wedge, about the mask
        % centroid; a wedge that runs off the frame does not sit at its bin centre
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
        % one wedge's contiguous run of populated bins, and its per-bin counts. Both
        % gate_wedge and accumulate call this, so verdict and block describe the same bins
            nB = obj.n_bin;
            sel = obj.tri.wedge == w & isfinite(obj.tri.bin) & isfinite(obj.tri.divergence);
            n_here = zeros(nB, 1);
            run_first = NaN;
            run_last = NaN;
            if ~any(sel)
                return
            end
            n_here = accumarray(obj.tri.bin(sel), 1, [nB 1]);
            % the run is over bins holding ANY triangle; min_tri_bin gates the means,
            % not the flux. Only an empty bin is a gap
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
