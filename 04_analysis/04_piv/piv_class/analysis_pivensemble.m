classdef analysis_pivensemble < handle
%ANALYSIS_PIVENSEMBLE  preprocess - interleave - correlate - gate, on the frames handed in
%   The caller cuts the frames (piv_getframescope); settings are param, set from
%   outside; the methods take only which result to compute.
%
% IN   stack_span  H x W x n_span numeric  from-halfwin : to+halfwin, before preprocessing
%      halfwin     1 x 1 int               frames either side of each endpoint
%      fps         1 x 1 double            frames per second
%      pixel2um    1 x 1 double            microns per pixel
%      clahe / clahe_size / wiener / wiener_size / prctile_low / prctile_high   piv_preprocess
% OUT  obj         handle                  endpoint and consecutive, one result struct each; see new_result
%
%   endpoint     lead n_frameset against trailing n_frameset : every pair spans the whole event
%   consecutive  pairs (1,2)(3,4).. one interval each, scaled by the intervals between the endpoints
%   USE          obj.correlate("endpoint")    obj.gate("endpoint")    maps = obj.get_corrplane("endpoint")

    properties
        % pivlab_ prefix : goes straight to PIVlab, theirs to define
        param = struct( ...
            'pivlab_windows',      [40 20; 20 10; 12 6], ... % P x 2 int  [win step] per pass
            'pivlab_mask_auto',    0, ...          % int         1 = suppress auto-corr peak
            'pivlab_subpixfinder', 1, ...          % int         1 = 3-point gauss, 2 = 2-D fit
            'pivlab_imdeform',     '*spline', ...  % char        deformation kernel, pass 2 on
            'pivlab_do_pad',       1, ...          % int         zero-pad -> linear not circular
            'pivlab_repeat',       1, ...          % int         pass 1 re-run on 4 quarter shifts
            'pivlab_neigh_thresh', [], ...         % float       gate_neighbour; [] = off
            'pivlab_mask',         [], ...         % H x W bool  true = pixel excluded
            'commonmode_tile',     40, ...         % float px    global-shift tile; [] = leave the shift in
            'tomasi_medfrac',      0.6, ...        % float       gate_tomasi;     [] = off
            'corr_floor_medfrac',  0.4, ...        % float       gate_corr_floor; [] = off
            'corr_floor_radius',   2, ...          % float px    plane cells this close to the peak are the peak
            'verbose',             true)           % bool        one line per gate, in -> out
    end

    properties (SetAccess = private)
        % settled at construction; changing one afterwards could not rebuild the stacks
        input = struct( ...
            'fps',      [], ...   % float      frames per second
            'pixel2um', [], ...   % float      microns per pixel
            'halfwin',  [], ...   % int        frames either side of each endpoint
            'n_span',   [], ...   % int        frames handed in
            'preproc',  [])       % struct     the pre-processing that was applied

        endpoint    % struct  total displacement across the event; see new_result
        consecutive % struct  the same total, walked one frame interval at a time
    end

    methods
        function obj = analysis_pivensemble(stack_span, halfwin, fps, pixel2um, opt)
            arguments
                stack_span {mustBeNumeric, mustBeNonempty}
                halfwin    (1,1) double {mustBeInteger, mustBeNonnegative}
                fps        (1,1) double {mustBePositive}
                pixel2um   (1,1) double {mustBePositive}
                opt.clahe        (1,1) logical = true   % lifts the sparse texture
                opt.clahe_size   (1,1) double  = 64     % tile size, px
                opt.wiener       (1,1) logical = true   % denoise; also the low-pass
                opt.wiener_size  (1,1) double  = 3      % px; sigma is half of it
                opt.prctile_low  (1,1) double  = 0.05   % imadjust, on the 0~1 frame
                opt.prctile_high (1,1) double  = 0.95
            end
            % 0. overlap check : the two framesets must not share a frame. The margin
            %    against the recording was the caller's (piv_getframescope)
            n_span     = size(stack_span, 3);
            n_frameset = 2*halfwin + 1;
            if n_span < 2*n_frameset
                error('analysis_pivensemble:framesetOverlap', ...
                    '%d frames with halfwin %d; the two framesets need %d', n_span, halfwin, 2*n_frameset);
            end
            % 1. the record
            obj.input.fps      = fps;
            obj.input.pixel2um = pixel2um;
            obj.input.halfwin  = halfwin;
            obj.input.n_span   = n_span;
            obj.input.preproc  = struct( ...
                'clahe',       opt.clahe,       'clahe_size',   opt.clahe_size, ...
                'wiener',      opt.wiener,      'wiener_size',  opt.wiener_size, ...
                'prctile_low', opt.prctile_low, 'prctile_high', opt.prctile_high);
            preproc_args = namedargs2cell(obj.input.preproc);

            % 2. preprocess once, from the record
            filtered = piv_preprocess(stack_span, preproc_args{:});

            % 3. consecutive : pairs (1,2)(3,4).. across the span, scaled by the intervals
            %    between the two endpoints
            n_pairs = floor(n_span/2);
            pair_consecutive = [(1:2:2*n_pairs)', (2:2:2*n_pairs)'];
            obj.consecutive = analysis_pivensemble.new_result(filtered, pair_consecutive, n_span - 1 - 2*halfwin);

            % 4. endpoint : (lead_1, trail_1, lead_2, trail_2, ...), every pair spans the event
            [interleaved, pair_endpoint] = piv_interleave(filtered, halfwin + 1, n_span - halfwin, halfwin);
            obj.endpoint = analysis_pivensemble.new_result(interleaved, pair_endpoint, 1);
        end

        function correlate(obj, mode)
        %CORRELATE  piv_corr_ensemble on a result's stack; the planes are kept without their maps
        %
        % IN   mode  "endpoint", "consecutive", or both (default)
        % OUT  none. result.planes gains xtable / ytable / utable / vtable / corr / typevector /
        %      imsize and peak_floor (ny x nx float, peak over the rms away from it, read off the
        %      final-pass planes before they are dropped -- see get_corrplane). xyuv stays []
        %      until gate() runs
            arguments
                obj
                mode (1,:) string {mustBeMember(mode, ["endpoint", "consecutive"])} = ["endpoint", "consecutive"]
            end
            for m = mode
                result = obj.(m);
                planes = obj.corr_planes(result.stack);
                planes.peak_floor = analysis_pivensemble.peak_floor_of(planes.maps, planes.typevector, obj.param.corr_floor_radius);
                planes.maps  = [];
                planes.corr2 = [];
                result.planes = planes;
                obj.(m) = result;
            end
        end

        function gate(obj, mode)
        %GATE  common mode out, then the gates in their order; every verdict kept, the kept set stamped
        %
        % IN   mode  "endpoint", "consecutive", or both (default). A mode not yet correlated is skipped
        % OUT  none. result gains xyuv, xyuv_ungated, uv, uv_ungated, common_mode, gates, gate_name
        %      gates  1 x 3 struct, one ROW per gate : name, stage (1 1 2; NaN = off in param),
        %             nVector_in, nVector_out (what it ALONE would leave), mask, info
        %
        %   SEQUENCE  common mode FIRST (a correction, not a gate). tomasi and corr_floor
        %   next, both stage 1 : one reads the images, the other the plane, neither the
        %   field. neighbour LAST : it reads its neighbours, so they must be clean already
            arguments
                obj
                mode (1,:) string {mustBeMember(mode, ["endpoint", "consecutive"])} = ["endpoint", "consecutive"]
            end
            for m = mode
                if isempty(obj.(m).planes)
                    if obj.param.verbose
                        fprintf('gate          skipping %s, not correlated\n', m);
                    end
                    continue
                end
                obj.(m) = obj.gate_result(obj.(m));
            end
        end

        function maps = get_corrplane(obj, mode)
        %GET_CORRPLANE  The final-pass correlation planes of one result, computed again and not kept
        %
        % IN   mode  "endpoint" | "consecutive"
        % OUT  maps  IA x IA x ny x nx x nPair single; sum(maps, 5) is the plane the peak was fitted to.
        %            Pair p's plane sits at (:, :, r, c, p) for the grid point planes.xtable(r, c)
            arguments
                obj
                mode (1,1) string {mustBeMember(mode, ["endpoint", "consecutive"])}
            end
            planes = obj.corr_planes(obj.(mode).stack);
            maps = planes.maps;
        end

        function [plane, info] = plane_at(obj, mode, maps, x, y, t)
        %PLANE_AT  One window's plane out of maps, summed over the pairs asked for
        %
        % IN   mode  "endpoint" | "consecutive", the result maps came from
        %      maps  from get_corrplane(mode)
        %      x, y  float  image coordinates, px; the nearest grid point is used
        %      t     int    which pairs to sum; omitted = every pair
        % OUT  plane IA x IA float  centre cell = zero displacement, the offset is the final-pass residual
        %      info  struct         mode, pairs, pair_frames, row, col, x, y, uv, valid
            arguments
                obj
                mode (1,1) string {mustBeMember(mode, ["endpoint", "consecutive"])}
                maps {mustBeNumeric}
                x (1,1) double
                y (1,1) double
                t double {mustBePositive, mustBeInteger} = []
            end
            result = obj.(mode);
            if isempty(result.planes)
                error('analysis_pivensemble:noPlanes', '%s has not been correlated', mode);
            end
            planes = result.planes;
            if isempty(t)
                t = 1:size(maps, 5);
            end
            dist = hypot(planes.xtable - x, planes.ytable - y);
            [~, nearest] = min(dist(:));
            [row, col]   = ind2sub(size(dist), nearest);
            plane = double(sum(maps(:, :, row, col, t), 5));
            if nargout > 1
                pix_row = round(planes.ytable(row, col));
                pix_col = round(planes.xtable(row, col));
                uv = [NaN NaN];
                if ~isempty(result.uv)
                    uv = squeeze(result.uv(pix_row, pix_col, :))';
                end
                info = struct('mode', mode, 'pairs', t, ...
                    'pair_frames', result.pair_frames(t, :), ...
                    'row', row, 'col', col, ...
                    'x', planes.xtable(row, col), 'y', planes.ytable(row, col), ...
                    'uv', uv, ...
                    'valid', planes.typevector(row, col) == 1);
            end
        end
    end

    methods (Static, Access = private)
        function result = new_result(stack, pair_frames, scale)
        % OUT  result  every field a result will ever carry, declared in one place
            result = struct( ...
                'stack',       stack, ...        % H x W x N float  filtered frames
                'pair_frames', pair_frames, ...  % P x 2 int        indices into the span
                'scale',       scale, ...        % int              onto uv -> event total
                'planes',      [], ...           % struct           piv_corr_ensemble's, maps dropped
                'xyuv',        [], ...           % ny x nx x 4 float  (:,:,1:2) [x y] window centre px,
                                ...              %   (:,:,3:4) [u v] UNSCALED, common mode out, GATED.
                                ...              %   The record; uv is the dense view of it
                'xyuv_ungated',[], ...           % ny x nx x 4 float  the same windows before the gates
                'common_mode', [0 0], ...        % 1 x 2 float      shift removed, scaled
                'gates',       struct([]), ...   % 1 x 3 struct     table; see gate()
                'gate_name',   "", ...           % string           the gates that were on, joined by +
                'uv_ungated',  [], ...           % H x W x 2 float  entering the gates
                'uv',          []);              % H x W x 2 float  leaving them, scaled
        end

        function ratio = peak_floor_of(maps, typevector, radius)
        %PEAK_FLOOR_OF  peak over the rms of the plane away from it, one number per window
        %   one peak on a quiet floor = structure; rival bumps = none. Final pass only,
        %   the one plane that is 1:1 with a vector
        %
        % IN   maps        IA x IA x ny x nx x nPair
        %      typevector  ny x nx        1 = a window that was correlated
        %      radius      1 x 1 float px cells this close to the peak are the peak
        % OUT  ratio       ny x nx float  NaN where nothing was correlated
            [ny, nx] = size(typevector);
            n = size(maps, 1);
            [YY, XX] = ndgrid(1:n, 1:n);
            ratio = NaN(ny, nx);
            for r = 1:ny
                for c = 1:nx
                    if typevector(r, c) ~= 1
                        continue
                    end
                    E = double(sum(maps(:, :, r, c, :), 5));
                    E = E - min(E(:));
                    [peak, k] = max(E(:));
                    if peak <= 0
                        continue
                    end
                    [pr, pc] = ind2sub(size(E), k);
                    away = hypot(XX - pc, YY - pr) > radius;
                    ratio(r, c) = peak / max(rms(E(away)), eps);
                end
            end
        end
    end

    methods (Access = private)
        function planes = corr_planes(obj, stack)
        % piv_corr_ensemble with this object's settings, maps included; correlate() and
        % get_corrplane() both come through here
            planes = piv_corr_ensemble(stack, ...
                window_sizes  = obj.param.pivlab_windows, ...
                exclmask      = obj.param.pivlab_mask, ...
                subpixfinder  = obj.param.pivlab_subpixfinder, ...
                mask_auto     = obj.param.pivlab_mask_auto, ...
                imdeform      = obj.param.pivlab_imdeform, ...
                repeat        = obj.param.pivlab_repeat, ...
                do_pad        = obj.param.pivlab_do_pad, ...
                save_corrmaps = true);
        end

        function result = gate_result(obj, result)
        % gate() for one result. See gate() for the sequence
            planes = result.planes;
            % 1. from the correlation, unscaled : the gates work in correlation px
            u = planes.utable;
            v = planes.vtable;

            % 2. global shift out. A correction, not a gate
            [u, v, common_mode] = obj.remove_commonmode(planes, u, v);
            result.common_mode  = common_mode * result.scale;
            % 2.1 one field with its own coordinates; the gates answer with masks, and the
            %     masks are applied at the two stamps
            xyuv     = cat(3, planes.xtable, planes.ytable, u, v);
            has_uv   = ~isnan(u) & ~isnan(v);
            keep_raw = (planes.typevector == 1) & has_uv;
            result.xyuv_ungated = piv_blank(xyuv, keep_raw);
            result.uv_ungated   = piv_stamp(xyuv, planes.imsize, keep_raw) * result.scale;

            % 3. stage 1 : two verdicts on the SAME field. A gate that is off in param
            %    answers with an all-false mask and info.on = false
            live_1 = has_uv;
            [mask_tomasi, info_tomasi] = obj.gate_tomasi(result, u, v);
            [mask_floor,  info_floor]  = obj.gate_corr_floor(result, u, v);
            rejected_1 = mask_tomasi | mask_floor;

            % 4. stage 2 : the only gate that reads the vectors around it, on a COPY with
            %    stage 1 applied
            u_stage2 = u;
            v_stage2 = v;
            u_stage2(rejected_1) = NaN;
            v_stage2(rejected_1) = NaN;
            live_2 = ~isnan(u_stage2) & ~isnan(v_stage2);
            [mask_neighbour, info_neighbour] = obj.gate_neighbour(result, u_stage2, v_stage2);
            rejected = rejected_1 | mask_neighbour;

            % 5. the table. stage NaN = off in param
            on = [info_tomasi.on, info_floor.on, info_neighbour.on];
            stage = [1 1 2];
            stage(~on) = NaN;
            n_1 = nnz(live_1);
            n_2 = nnz(live_2);
            result.gates = struct( ...
                'name',        {"tomasi", "corr_floor", "neighbour"}, ...
                'stage',       num2cell(stage), ...
                'nVector_in',  {n_1, n_1, n_2}, ...
                'nVector_out', {n_1 - nnz(mask_tomasi    & live_1), ...
                                n_1 - nnz(mask_floor     & live_1), ...
                                n_2 - nnz(mask_neighbour & live_2)}, ...
                'mask',        {mask_tomasi, mask_floor, mask_neighbour}, ...
                'info',        {info_tomasi, info_floor, info_neighbour});
            gate_list = ["tomasi", "corr_floor", "neighbour"];
            gate_list = gate_list(on);
            if ~isempty(obj.param.commonmode_tile)
                gate_list = ["commonmode", gate_list];
            end
            result.gate_name = strjoin(gate_list, "+");

            % 6. onto the image, scaled to the event total. The same xyuv as 2.1, a different mask
            keep_gated   = keep_raw & ~rejected;
            result.xyuv  = piv_blank(xyuv, keep_gated);
            result.uv    = piv_stamp(xyuv, planes.imsize, keep_gated) * result.scale;

            if obj.param.verbose
                for g = result.gates
                    note = "";
                    if ~g.info.on
                        note = "off";
                    end
                    fprintf('  %-3g %-10s %4d -> %-4d  %s\n', g.stage, g.name, g.nVector_in, g.nVector_out, note);
                end
                len = hypot(result.uv(:,:,1), result.uv(:,:,2));
                fprintf('%15s %4d left | |d| p50 %.3f px = %.3f um | %s\n', '', ...
                    nnz(~isnan(len)), median(len(:), 'omitnan'), ...
                    median(len(:), 'omitnan') * obj.input.pixel2um, result.gate_name);
            end
        end

        function [u, v, common_mode] = remove_commonmode(obj, planes, u, v)
        % IN   planes  struct         for the grid spacing : tile read in px
        % OUT  u, v    ny x nx float  the field with the global shift taken out
        %      common_mode  1 x 2 float  what was removed, UNSCALED. [0 0] when the tile is empty
        %
        %   tile medians, then their median : weighted by area. The amount is returned,
        %   not discarded -- drift and tissue translation are not separable here
            common_mode = [0 0];
            if isempty(obj.param.commonmode_tile)
                return
            end
            grid_step  = planes.xtable(1, 2) - planes.xtable(1, 1);
            has_vector = ~isnan(u) & ~isnan(v);
            % 1. tile index of every grid point
            tile_cells = max(1, round(obj.param.commonmode_tile / grid_step));
            [ny, nx]   = size(u);
            tile_row   = ceil((1:ny)'/tile_cells) * ones(1, nx);
            tile_col   = ones(ny, 1) * ceil((1:nx)/tile_cells);
            tile_id    = (tile_row - 1)*max(tile_col(:)) + tile_col;
            % 2. one median vector per tile; too-sparse tiles drop out
            n_in_tile    = accumarray(tile_id(has_vector), 1);
            tile_u       = accumarray(tile_id(has_vector), u(has_vector), [], @median);
            tile_v       = accumarray(tile_id(has_vector), v(has_vector), [], @median);
            dense_enough = n_in_tile >= 20;      % never met at grid step 10 px; see CLAUDE_LOG.md
            if ~any(dense_enough)
                return
            end
            common_mode = [median(tile_u(dense_enough)), median(tile_v(dense_enough))];
            u = u - common_mode(1);
            v = v - common_mode(2);
        end

        function [mask, info] = gate_tomasi(obj, result, u, ~)
        % IN   result  struct        read for the IMAGES and the grid, never the field
        %      u       ny x nx float for the grid SIZE only
        % OUT  mask    ny x nx bool  true = nothing here to track
        %      info    struct        on, tomasi_medfrac, median, cut,
        %                            lambda2 (ny x nx x 2 float, (:,:,1) = leading end)
        %
        %   from the images, not the correlation. Both ends, judged on the weaker; the
        %   cut is a fraction of this frame's own median
            info = struct('on', ~isempty(obj.param.tomasi_medfrac), ...
                'tomasi_medfrac', obj.param.tomasi_medfrac, ...
                'median', NaN, 'cut', NaN, 'lambda2', []);
            mask = false(size(u));
            if ~info.on
                return
            end
            planes    = result.planes;
            window_px = obj.param.pivlab_windows(end, 1);
            lead_img  = mean(result.stack(:, :, 1:2:end), 3);
            trail_img = mean(result.stack(:, :, 2:2:end), 3);
            info.lambda2 = cat(3, ...
                piv_trackability(lead_img,  planes.xtable, planes.ytable, window_px, obj.param.pivlab_mask), ...
                piv_trackability(trail_img, planes.xtable, planes.ytable, window_px, obj.param.pivlab_mask));
            worst_end = min(info.lambda2, [], 3);
            info.median = median(worst_end(:), 'omitnan');
            info.cut    = obj.param.tomasi_medfrac * info.median;
            % NaN means piv_trackability could not build a tensor : the window was almost
            % entirely masked out. Nothing to track, so out
            mask = worst_end <= info.cut | isnan(worst_end);
        end

        function [mask, info] = gate_corr_floor(obj, result, u, ~)
        % IN   result  read for planes.peak_floor, never for the field
        %      u       for the grid size
        % OUT  mask    ny x nx bool  true = the plane has no clean peak
        %      info    struct        on, corr_floor_medfrac, median, cut, peak_floor (ny x nx)
        %
        %   the statistic comes from correlate(), so no maps are needed here
            info = struct('on', ~isempty(obj.param.corr_floor_medfrac), ...
                'corr_floor_medfrac', obj.param.corr_floor_medfrac, ...
                'median', NaN, 'cut', NaN, 'peak_floor', []);
            mask = false(size(u));
            if ~info.on
                return
            end
            ratio = result.planes.peak_floor;
            info.peak_floor = ratio;
            info.median     = median(ratio(:), 'omitnan');
            info.cut        = obj.param.corr_floor_medfrac * info.median;
            mask = ratio <= info.cut | isnan(ratio);
        end

        function [mask, info] = gate_neighbour(obj, ~, u, v)
        % OUT  mask  ny x nx bool  true = vector disagrees with its neighbourhood
        %      info  struct        on, pivlab_neigh_thresh
        %
        %   Westerweel-Scarano universal median. Reads the FIELD, hence stage 2
            info = struct('on', ~isempty(obj.param.pivlab_neigh_thresh), ...
                'pivlab_neigh_thresh', obj.param.pivlab_neigh_thresh);
            mask = false(size(u));
            if ~info.on
                return
            end
            [filt_u, filt_v] = postproc.PIVlab_postproc(u, v, 1, 1, [], false, ...
                [], true, obj.param.pivlab_neigh_thresh);
            mask = (isnan(filt_u) | isnan(filt_v)) & ~(isnan(u) | isnan(v));
        end
    end
end
