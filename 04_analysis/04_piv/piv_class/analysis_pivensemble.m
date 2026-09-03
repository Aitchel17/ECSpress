classdef analysis_pivensemble < handle
% Ensemble-correlation PIV across the two endpoints of one dilation event.
%
% IN   stack_span H x W x n_span   from-halfwin : to+halfwin, before preprocessing;
%                                  piv_getframescope cuts it out of the recording
%      halfwin    int              frames either side of each endpoint
%      fps        float            frames per second
%      pixel2um   float            microns per pixel
%      name-value                  pre-processing settings, applied once here
% OUT  obj        handle           knows the frames it was handed. Nothing about the
%                                  recording, events, sleep state or the session
%
% RESULTS  both come out as TOTAL px across the event
%   endpoint     lead 2*halfwin+1 vs trailing : every pair spans the whole event
%   consecutive  pairs (1,2)(3,4).. one interval each : x interval count
%   why both     consecutive reaching endpoint's answer : it accumulated
%   err          a jump : all but one plane at zero, plane sum votes it to nothing
%   caution      scale is applied inside : nothing outside ever multiplies a field
%
% USE
%   obj.correlate()             both, ~11 s. PIVlab multipass, irreversible
%   obj.correlate("endpoint")   one only, most of the cost saved
%   obj.gate()                  ~0.09 s. Always restarts from the correlation
%   obj.gate(tomasi = false)    the field without one gate, threshold kept
%   obj.plane_at(x, y[, t])     the correlation plane behind one vector
%
% DESIGN
%   split        correlate expensive, gate cheap : retry a gate for free
%   gates        each returns [mask, info] : judges, changes nothing
%   order        is a finding : gate() spells the sequence out, never loops it
%   no facade    postproc.PIVlab_* direct, not piv_validate : a facade bundling
%                three stages takes the ordering away from the class that owns it
%   margin       the caller cuts from-halfwin : to+halfwin (piv_getframescope);
%                only the overlap check lives here
%   blink        endpoint.stack interleaved : sliceViewer flips FROM/TO

    properties
        % pivlab_ prefix : goes straight to PIVlab, theirs to define
        % no prefix      : defined here
        % all read by gate() at call time : change one, re-gate, no re-correlation
        %
        % THREE GATES, on 15 hand-labelled window groups across two events of
        % opposite polarity (ev2 dilation, ev22 constriction). Every setting below
        % keeps 100% of the labelled-good windows.
        %
        % tomasi_medfrac -- Shi-Tomasi lambda2 as a fraction of the frame's median
        %   why       the failure is in the IMAGE, not the match : a window with no
        %             corner still correlates sharply, and its peak then sits
        %             anywhere along the flat direction
        %   why       relative, not absolute : lambda2's median moved 1.25x between
        %             two events of the same recording, and a fixed 1.6e-4 that cut
        %             21% on one cut only 8% on the other, letting labelled-bad
        %             windows through
        %   why       relative, not a percentile : a percentile removes its share
        %             even on a frame where nothing is wrong. Measured, the SHAPE
        %             of the distribution matches between events (p20/p50 = 0.593
        %             vs 0.642) while the scale does not, which is exactly what a
        %             ratio-to-median follows
        %   see FINDINGS.md
        %   caution   the median holds only while bad windows are a MINORITY
        %
        % corr_floor_medfrac -- peak height over the rms of the plane AWAY from it
        %   why       a window with real structure gives one peak on a clean floor;
        %             one without gives a floor full of rival bumps. Read on the
        %             final pass, the only plane that is 1:1 with a vector
        %   see FINDINGS.md
        %   caution   weaker than lambda2 at every budget. 0.4 cuts 13% and keeps
        %             100% of labelled-good; 0.5 already costs 2-9%
        %
        % pivlab_neigh_thresh -- Westerweel-Scarano. OFF by default, see gate()
        %   see FINDINGS.md
        %   why       it is the only gate that reads a vector's NEIGHBOURS, so a
        %             coherent patch moving faster than the field around it is
        %             exactly what it is built to remove -- and that patch is the
        %             signal here
        %   note      the threshold stays at PIVlab's 2 : the switch is off, not the
        %             gate, so gates(3).mask still shows what it would have taken.
        %             gate(neighbour = true) puts it back for one call
        %
        % RETIRED, in zz_notinuse/piv_legacy/ : gate_rheight (corr2 reads displacement, not
        % failure, corr(|d|,corr2) = -0.475) and gate_rwidth (peak width in an
        % ensemble is displacement spread, and it had a direction preference)
        %
        % THE WINDOW LADDER, settled 260803 on HQL99 260601_005. Moved here from
        % piv_integration_testbed's header when that file was replaced 260808 --
        % these justify the CORRELATION SETTINGS, so they belong beside the
        % setting rather than in a script that happened to measure them.
        %
        % pivlab_windows -- [window step] PER PASS, one row each
        %   see FINDINGS.md
        %   err       more ENSEMBLE PAIRS does not fix a wide peak : 5 -> 53 pairs
        %             moved the half-width 1.99 -> 1.93 px only. The width is window
        %             geometry, not frame noise, so it is bought here and nowhere
        %             else. The event span caps the pair count anyway -- framesets
        %             must not overlap
        %   caution   corr2 cannot see this failure. Blur holds it at 0.73 at the
        %             wall while the peak spreads underneath, which is what the
        %             corr_floor gate above is for
        %
        % coremask is NOT grown before the polar measures read from it
        %   why       the traced boundary IS the ring origin : it was drawn on the
        %             maximum-dilation frame, so growing it pushes ring 1 into
        %             tissue the PVS never reaches
        %   see FINDINGS.md
        param = struct( ...
            'pivlab_windows',      [40 20; 20 10; 12 6], ... % P x 2 int  [win step] per pass
            'pivlab_mask_auto',    0, ...          % int         1 = suppress auto-corr peak
            'pivlab_subpixfinder', 1, ...          % int         1 = 3-point gauss, 2 = 2-D fit
            'pivlab_imdeform',     '*spline', ...  % char        deformation kernel, pass 2 on
            'pivlab_do_pad',       1, ...          % int         zero-pad -> linear not circular
            'pivlab_repeat',       1, ...          % int         pass 1 re-run on 4 quarter shifts
            'pivlab_neigh_thresh', 2, ...          % float       gate_neighbour; [] = off
            'pivlab_mask',         [], ...         % H x W bool  true = pixel excluded
            'commonmode_tile',     40, ...         % float px    global-shift tile; [] leaves
                                     ...           %             the shift in, which lets a
                                     ...           %             gate pick surviving sectors
            'tomasi_medfrac',      0.6, ...        % float       gate_tomasi;     [] = off
            'corr_floor_medfrac',  0.4, ...        % float       gate_corr_floor; [] = off
            'verbose',             true)           % bool        one line per gate, in -> out.
                                                   %             No report method: a row that
                                                   %             looks wrong -> gates(k).info
    end

    properties (SetAccess = private)
        % settled at construction : each one decides which frames get prepared, so
        % changing them afterwards could not rebuild the stacks
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
        % why  pre-processing is contrast only, nothing that moves a pixel
        % why  CLAHE lifts the sparse texture the correlation keys on, Wiener
        %      takes the shot noise with it
        % why  here and not a method : it has to happen exactly once, first
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
            % 0. The class's own contract : the two framesets must not overlap, or a
            %    frame correlates with itself. Where the span sits in the recording is
            %    not checked here -- the caller cut it
            n_span     = size(stack_span, 3);
            n_frameset = 2*halfwin + 1;      % frames making up one endpoint's set
            if n_span < 2*n_frameset
                error('analysis_pivensemble:framesetOverlap', ...
                    '%d frames with halfwin %d; the two framesets need %d', ...
                    n_span, halfwin, 2*n_frameset);
            end

            obj.input.fps      = fps;
            obj.input.pixel2um = pixel2um;
            obj.input.halfwin  = halfwin;
            obj.input.n_span   = n_span;

            % 1. Pre-processing, listed POSITIVELY rather than subtracted
            %    why      halfwin is not pre-processing : it picks which frames
            %             exist, not how they are filtered. It has input.halfwin
            %    caution  rmfield(opt,...) would let a future non-filter option
            %             land in a record that claims to say what was applied
            obj.input.preproc = struct( ...
                'clahe',       opt.clahe,       'clahe_size',   opt.clahe_size, ...
                'wiener',      opt.wiener,      'wiener_size',  opt.wiener_size, ...
                'prctile_low', opt.prctile_low, 'prctile_high', opt.prctile_high);

            % 2. Filter once, FROM that record
            %    why  record is the source of the call, not a copy : cannot drift
            %    why  raw stack not kept : the two schemes cannot diverge
            preproc_args = namedargs2cell(obj.input.preproc);
            filtered     = piv_preprocess(stack_span, preproc_args{:});

            % 3. Consecutive : the span as it stands, pairs 1:2:end
            %    caution  scaled by INTERVALS between the endpoints (to - from), not
            %             by pair count : the pairs sample every other interval, so
            %             counting pairs would halve it
            n_pairs = floor(n_span/2);
            pair_consecutive = [(1:2:2*n_pairs)', (2:2:2*n_pairs)'];
            obj.consecutive = analysis_pivensemble.new_result(filtered, ...
                pair_consecutive, n_span - 1 - 2*halfwin);

            % 4. Endpoint : (lead_1, trail_1, lead_2, trail_2, ...)
            %    why  every pair already spans the event, so scale 1. pair_frames
            %         are indices into the span; the caller's span maps them back
            [interleaved, pair_endpoint] = piv_interleave(filtered, ...
                halfwin + 1, n_span - halfwin, halfwin);
            obj.endpoint = analysis_pivensemble.new_result(interleaved, pair_endpoint, 1);
        end

        function correlate(obj, which, opt)
        % IN   which         str     "endpoint", "consecutive", or both (default)
        %      window_sizes  P x 2 int  one ROW PER PASS, coarse to fine, each row
        %                    [window step]. [40 20; 20 10] is two passes,
        %                    [40 20; 20 10; 12 6] is three. Omitted = obj.param
        %                    caution  only row 1's STEP is used. piv_corr_ensemble
        %                    passes one step_size to PIVlab and derives the later
        %                    passes from the interrogation areas, so writing
        %                    [40 20; 20 5] does not give pass 2 a step of 5
        % OUT  none; fills planes. xyuv stays [] until gate() runs -- it means the
        %      GATED grid and nothing else, so there is no window in which it
        %      means the raw one
        %
        % why      the expensive half and the only irreversible one
        % why      gate() separate : it is the half worth repeating, and doing it
        %          here would gate a result the caller may not have correlated
        % why      piv_corr_ensemble gates nothing : every judgement is in gate(),
        %          in an order this class chose
        % why      planes always kept : that is what makes gate() cheap
        % why      a ladder given here is WRITTEN INTO obj.param before it runs.
        %          The planes are the only record of what was correlated, so an
        %          argument that left the property saying something else would
        %          make the object disagree with its own contents
        % caution  name one result to skip the other, which is most of the cost
            arguments
                obj
                which (1,:) string {mustBeMember(which, ...
                    ["endpoint", "consecutive"])} = ["endpoint", "consecutive"]
                opt.window_sizes (:,2) double {mustBeInteger, mustBePositive} = zeros(0,2)
            end
            if ~isempty(opt.window_sizes)
                obj.param.pivlab_windows = opt.window_sizes;
            end
            settings = obj.param;
            for name = which
                result = obj.(name);
                planes = piv_corr_ensemble(result.stack, ...
                    window_sizes  = settings.pivlab_windows, ...
                    exclmask      = settings.pivlab_mask, ...
                    subpixfinder  = settings.pivlab_subpixfinder, ...
                    mask_auto     = settings.pivlab_mask_auto, ...
                    imdeform      = settings.pivlab_imdeform, ...
                    repeat        = settings.pivlab_repeat, ...
                    do_pad        = settings.pivlab_do_pad, ...
                    save_corrmaps = true);
                result.planes = planes;
                % caution  UNSCALED. result.scale goes on in gate() step 6 only
                % why      gate thresholds are in correlation px : the retired gate_rwidth
                %          measures a vector against its own peak half-width, and
                %          a consecutive field already x54 clears that bar for
                %          nothing
                obj.(name)     = result;
            end
        end

        function gate(obj, which, opt)
        % IN   which   "endpoint", "consecutive", or both (default)
        %      tomasi / corr_floor   bool, default TRUE
        %      neighbour             bool, default FALSE -- it rejects labelled-good
        %              vectors 3x more often than labelled-bad ones; see param
        %              a switched-off gate still forms its verdict into gates(k).mask
        % OUT  none; fills xyuv, xyuv_ungated, uv, uv_ungated, common_mode, gates
        %
        % why      switch vs param = [] : the switch keeps the THRESHOLD, so the
        %          field without a gate can be seen without losing the number
        % caution  restarts from the correlation every call : 2e-4 then 1.6e-4
        %          gives exactly what 1.6e-4 alone gives, calls never stack
        %
        % result.gates -- 1 x 4 struct, one ROW per gate:
        %   name         str           which gate
        %   stage        int           pass it ran in, 1 1 1 2. NaN = switch was off
        %   nVector_in   int           vectors alive when it formed its verdict
        %   nVector_out  int           what it ALONE would leave
        %   mask         ny x nx bool  its verdict, whether carried out or not
        %   info         struct        thresholds given + what it measured
        % why      stage is the dependency as a VALUE, not as the order the code
        %          happens to be written in
        % caution  in/out describe that gate ALONE : rows sharing a stage are NOT
        %          cumulative -- same field in, overlapping windows out
        % note     info.on false : the threshold in param was empty
        %
        % SEQUENCE  the order is the result, so it is written out
        %   common mode  FIRST. Left in, every vector is longer one way and
        %                shorter the other : a length-keyed gate turns directional
        %   tomasi       next. The only judgement read from the IMAGES, so it
        %                cannot be confounded by the correlation it is filtering
        %   corr_floor   same stage as tomasi : it reads the plane, not the vectors
        %   neighbour    LAST, when asked for. Reads its neighbours, so they must be
        %                clean already -- which is why it stays stage 2 even off
        % caution  the common mode is NOT switchable here : it is a correction,
        %          not a gate, and uv_ungated is stamped from it.
        %          param.commonmode_tile = [] turns it off
            arguments
                obj
                which (1,:) string {mustBeMember(which, ...
                    ["endpoint", "consecutive"])} = ["endpoint", "consecutive"]
                opt.tomasi     (1,1) logical = true
                opt.corr_floor (1,1) logical = true
                opt.neighbour  (1,1) logical = false
            end
            % Gate what has actually been correlated. correlate() takes a name to
            % skip the other half's cost, and making the caller repeat that name
            % here is a duplication that drifts -- correlate one and gate both is
            % the trap the pairing sets. Skipped is reported, never silent
            skipped = which(arrayfun(@(n) isempty(obj.(n).planes), which));
            which   = setdiff(which, skipped, 'stable');
            if isempty(which)
                error('analysis_pivensemble:noPlanes', ...
                    'nothing has been correlated; call correlate() first.');
            end
            if ~isempty(skipped) && obj.param.verbose
                fprintf('gate          skipping %s, not correlated\n', ...
                    strjoin(skipped, ', '));
            end

            for name = which
                result = obj.(name);
                planes = result.planes;
                % 1. From the correlation, unscaled : the gates work in
                %    correlation px, see step 6
                u = planes.utable;
                v = planes.vtable;

                % 2. Global shift out. A correction, not a gate
                %    see FINDINGS.md
                [u, v, common_mode] = obj.remove_commonmode(planes, u, v);
                result.common_mode  = common_mode * result.scale;
                % 2.1 One field, carried with its own coordinates. u and v are not
                %     written to again below -- every gate answers with a mask and
                %     the masks are applied at the two stamps
                xyuv     = cat(3, planes.xtable, planes.ytable, u, v);
                has_uv   = ~isnan(u) & ~isnan(v);
                keep_raw = (planes.typevector == 1) & has_uv;
                result.xyuv_ungated = piv_blank(xyuv, keep_raw);
                result.uv_ungated   = piv_stamp(xyuv, planes.imsize, keep_raw) * result.scale;
                if obj.param.verbose
                    fprintf('%-12s %4d raw | common mode (%+.3f, %+.3f) px out\n', ...
                        name, nnz(~isnan(planes.utable)), result.common_mode);
                end

                % 3. Stage 1 : two verdicts on the SAME field
                %    why  neither reads a vector's neighbours -- tomasi reads the
                %         IMAGES, corr_floor the plane around its own peak -- so
                %         their order cannot matter and their masks just stack
                %    why  every gate answers even when switched off : gates(k).mask
                %         then shows what it WOULD have taken
                live_1 = has_uv;
                [mask_tomasi, info_tomasi] = obj.gate_tomasi(result, u, v);
                [mask_floor,  info_floor]  = obj.gate_corr_floor(result, u, v);

                % 3.1 Collect the switched-on ones. A verdict, nothing applied yet
                rejected_1 = false(size(u));
                if opt.tomasi;     rejected_1 = rejected_1 | mask_tomasi; end
                if opt.corr_floor; rejected_1 = rejected_1 | mask_floor;  end

                % 4. Stage 2 : the only gate that reads the vectors around it
                %    err  on the uncleaned field, a vector the others rejected
                %         still votes in its neighbour's median
                %    why  stage 1 is applied to a COPY. u and v stay as the
                %         correlation left them, so uv_ungated and uv come off the
                %         same field and differ only in their mask
                u_stage2 = u;
                v_stage2 = v;
                u_stage2(rejected_1) = NaN;
                v_stage2(rejected_1) = NaN;
                live_2 = ~isnan(u_stage2) & ~isnan(v_stage2);
                [mask_neighbour, info_neighbour] = obj.gate_neighbour(result, u_stage2, v_stage2);
                rejected = rejected_1;
                if opt.neighbour
                    rejected = rejected | mask_neighbour;
                end

                % 5. The table. stage NaN = the switch was off
                %    see FINDINGS.md
                stage = [1 1 2];
                stage(~[opt.tomasi, opt.corr_floor, opt.neighbour]) = NaN;
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

                % 6. Onto the image, scaled to the event total
                %    caution  scale is applied HERE and only here
                %    why  the same xyuv as step 2.1, a different mask. The gated
                %         and ungated maps cannot disagree about which vector sat
                %         where, because neither was ever written to
                keep_gated   = keep_raw & ~rejected;
                result.xyuv  = piv_blank(xyuv, keep_gated);
                result.uv    = piv_stamp(xyuv, planes.imsize, keep_gated) * result.scale;
                obj.(name)   = result;

                if obj.param.verbose
                    for g = result.gates
                        note = "";
                        if ~g.info.on;         note = "off";
                        elseif isnan(g.stage); note = "skipped";
                        end
                        fprintf('  %-3g %-10s %4d -> %-4d  %s\n', g.stage, ...
                            g.name, g.nVector_in, g.nVector_out, note);
                    end
                    len = hypot(result.uv(:,:,1), result.uv(:,:,2));
                    fprintf('%15s %4d left | |d| p50 %.3f px = %.3f um\n', '', ...
                        nnz(~isnan(len)), median(len(:), 'omitnan'), ...
                        median(len(:), 'omitnan') * obj.input.pixel2um);
                end
            end
        end

        function [plane, info] = plane_at(obj, x, y, t)
        % IN   x, y  float  image coordinates in px; nearest grid point is used
        %      t     int    which pairs to sum, and by naming them, which result:
        %              omitted  ENDPOINT ensemble, every pair -- the plane the
        %                       reported displacement was fitted to
        %              7        CONSECUTIVE, one frame interval
        %              1:5      CONSECUTIVE, the first five summed
        % OUT  plane IA x IA float  centre cell = zero displacement, so the peak
        %                          offset is the final-pass residual
        %      info  struct         which result, which pairs, grid cell, the
        %                          coordinates used, the vector there, validity
        %
        % why  instance not static : nothing to look at until correlate() has run
        % why  planes captured before any gating : a window the gates threw away
        %      still has one, which is the point of looking
            arguments
                obj
                x (1,1) double
                y (1,1) double
                t double {mustBePositive, mustBeInteger} = []
            end
            if isempty(t); name = "endpoint"; else; name = "consecutive"; end
            result = obj.(name);
            if isempty(result.planes)
                error('analysis_pivensemble:noPlanes', ...
                    '%s has not been correlated; call correlate("%s").', name, name);
            end
            planes = result.planes;
            if isempty(t); t = 1:size(planes.maps, 5); end

            dist = hypot(planes.xtable - x, planes.ytable - y);
            [~, nearest] = min(dist(:));
            [row, col]   = ind2sub(size(dist), nearest);
            plane = double(sum(planes.maps(:, :, row, col, t), 5));
            if nargout > 1
                pix_row = round(planes.ytable(row, col));
                pix_col = round(planes.xtable(row, col));
                uv = [NaN NaN];
                if ~isempty(result.uv)
                    uv = squeeze(result.uv(pix_row, pix_col, :))';
                end
                info = struct('which', name, 'pairs', t, ...
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
        %
        % why  scale is structural : fixed by the pairing at construction, so no
        %      caller and no method decides it later
            result = struct( ...
                'stack',       stack, ...        % H x W x N float  filtered frames
                'pair_frames', pair_frames, ...  % P x 2 int        indices into the span
                'scale',       scale, ...        % int              onto uv -> event total
                'planes',      [], ...           % struct           piv_corr_ensemble's
                'xyuv',        [], ...           % ny x nx x 4 float  (:,:,1:2) [x y]
                                ...              %   window centre px, (:,:,3:4) [u v]
                                ...              %   UNSCALED, common mode out, GATED.
                                ...              %   The record; uv is the dense view
                                ...              %   of it. [] until gate() has run
                'xyuv_ungated',[], ...           % ny x nx x 4 float  the same windows
                                ...              %   before the gates. Same name and
                                ...              %   same meaning as piv_run_events'''s
                'common_mode', [0 0], ...        % 1 x 2 float      shift removed, scaled
                'gates',       struct([]), ...   % 1 x 4 struct     table; see gate()
                'uv_ungated',  [], ...           % H x W x 2 float  entering the gates
                'uv',          []);              % H x W x 2 float  leaving them, scaled
        % note  the dense corr map is gone. It was written here and in three other
        %       places and read in none; the corr that gets used is planes.corr,
        %       one value per window, and it is consumed before the stamp
        end
    end

    methods (Access = private)
        function [u, v, common_mode] = remove_commonmode(obj, planes, u, v)
        % IN   planes  struct         for the grid spacing : tile read in px
        % OUT  u, v    ny x nx float  the field with the global shift taken out
        %      common_mode  1 x 2 float  what was removed, UNSCALED. [0 0] when
        %                                the tile is empty
        %
        % why      the one thing that CHANGES a vector rather than removing it,
        %          so it is not a gate and has no switch in gate()
        % why      radially symmetric dilation cancels to zero : what is left over
        %          the whole field is a rigid shift
        % caution  registration drift vs real tissue translation : not separable
        %          from one event, so common_mode is handed back, not discarded
        % see FINDINGS.md
        % err      a plain median then leans toward wherever the field is
        %          strongest -- one median per tile, then median of tiles, weights
        %          by AREA instead
            common_mode = [0 0];
            if isempty(obj.param.commonmode_tile); return; end

            grid_step  = planes.xtable(1, 2) - planes.xtable(1, 1);
            has_vector = ~isnan(u) & ~isnan(v);
            % 1. Tile index of every grid point
            tile_cells = max(1, round(obj.param.commonmode_tile / grid_step));
            [ny, nx]   = size(u);
            tile_row   = ceil((1:ny)'/tile_cells) * ones(1, nx);
            tile_col   = ones(ny, 1) * ceil((1:nx)/tile_cells);
            tile_id    = (tile_row - 1)*max(tile_col(:)) + tile_col;
            % 2. One median vector per tile; too-sparse tiles drop out
            n_in_tile    = accumarray(tile_id(has_vector), 1);
            tile_u       = accumarray(tile_id(has_vector), u(has_vector), [], @median);
            tile_v       = accumarray(tile_id(has_vector), v(has_vector), [], @median);
            dense_enough = n_in_tile >= 20;
            if ~any(dense_enough); return; end

            common_mode = [median(tile_u(dense_enough)), median(tile_v(dense_enough))];
            u = u - common_mode(1);
            v = v - common_mode(2);
        end

        function [mask, info] = gate_tomasi(obj, result, u, ~)
        % IN   result  struct  read for the IMAGES and the grid, never the field
        %      u       ny x nx float  for the grid SIZE only : this gate never
        %                             looks at the vectors, which is the point
        % OUT  mask  ny x nx bool  true = nothing here to track
        %      info  struct        on(bool) tomasi_lambda2(float)
        %                          lambda2(ny x nx x 2 float, (:,:,1) = leading end)
        %
        % why      read from the images : owes nothing to the correlation it judges
        % why      BOTH ends, judged on the weaker : a window that had a corner and
        %          lost it by the second frame is no more trackable than one that
        %          never had one
        % note     lambda2(:,:,1) = leading frames, (:,:,2) = trailing, so which
        %          end failed stays visible afterwards
        % caution  a percentile : removes its share whether or not anything is bad.
        %          A choice about how much to spend, not a test anything failed
            info = struct('on', ~isempty(obj.param.tomasi_medfrac), ...
                'tomasi_medfrac', obj.param.tomasi_medfrac, ...
                'median', NaN, 'cut', NaN, 'lambda2', []);
            mask = false(size(u));
            if ~info.on; return; end

            planes    = result.planes;
            window_px = obj.param.pivlab_windows(end, 1);
            lead_img  = mean(result.stack(:, :, 1:2:end), 3);
            trail_img = mean(result.stack(:, :, 2:2:end), 3);
            info.lambda2 = cat(3, ...
                piv_trackability(lead_img,  planes.xtable, planes.ytable, ...
                    window_px, obj.param.pivlab_mask), ...
                piv_trackability(trail_img, planes.xtable, planes.ytable, ...
                    window_px, obj.param.pivlab_mask));
            worst_end = min(info.lambda2, [], 3);
            % The cut rides on this frame's own median, so it follows the scale
            % without forcing a fixed share out. median, not mean: it holds as
            % long as the bad windows are a MINORITY, which is the assumption
            info.median = median(worst_end(:), 'omitnan');
            info.cut    = obj.param.tomasi_medfrac * info.median;
            % NaN means piv_trackability could not build a tensor at all -- the
            % window was almost entirely masked out. Nothing to track, so out
            mask = worst_end <= info.cut | isnan(worst_end);
        end

        function [mask, info] = gate_corr_floor(obj, result, u, ~)
        % IN   result  read for the PLANES only, never for the field
        %      u       for the grid size; the vectors themselves are not used
        % OUT  mask    ny x nx bool  true = the plane has no clean peak
        %      info    struct        on(bool) corr_floor_medfrac(float)
        %                            median(float) cut(float)
        %                            peak_floor(ny x nx float, the ratio itself)
        %
        % why      a window with real structure gives ONE peak on a quiet floor;
        %          one without gives a floor full of rival bumps. The information
        %          is in what surrounds the peak, not in the peak
        % why      the final pass, because it is the only plane that is 1:1 with a
        %          vector. The earlier passes cover 3-12 vectors each
        % why      rms of everything more than 2 px from the peak : 2 px is where
        %          the peak itself has decayed, measured on this data
        % see FINDINGS.md
        %
        % DIRECTION -- the disc is isotropic, and that is NOT settled
        % note     the obvious refinement is to cut the peak's own ridge out before
        %          taking the rms. Tested on ev2, and the ridge is PERPENDICULAR to
        %          the displacement, not along it : a linear feature correlates into
        %          a ridge running along the FEATURE, and the feature is at right
        %          angles to the motion it produces. The aperture problem, which is
        %          also why gate_tomasi wins
        % see FINDINGS.md
        % caution  ONE EVENT. gate_tomasi was settled on ev2 AND ev22 and the
        %          absolute-threshold version looked just as good on ev2 alone.
        %          Not switched until ev22 says the same
        % caution  weaker than gate_tomasi at every budget, and nearly independent
        %          of it (22-29% overlap). It is a second opinion, not a
        %          replacement
            info = struct('on', ~isempty(obj.param.corr_floor_medfrac), ...
                'corr_floor_medfrac', obj.param.corr_floor_medfrac, ...
                'median', NaN, 'cut', NaN, 'peak_floor', []);
            mask = false(size(u));
            if ~info.on || isempty(result.planes.maps); return; end

            planes = result.planes;
            [ny, nx] = size(u);
            n = size(planes.maps, 1);
            [YY, XX] = ndgrid(1:n, 1:n);
            ratio = NaN(ny, nx);
            for r = 1:ny
                for c = 1:nx
                    if planes.typevector(r,c) ~= 1; continue; end
                    E = double(sum(planes.maps(:, :, r, c, :), 5));
                    E = E - min(E(:));
                    [peak, k] = max(E(:));
                    if peak <= 0; continue; end
                    [pr, pc] = ind2sub(size(E), k);
                    away = hypot(XX - pc, YY - pr) > 2;
                    ratio(r,c) = peak / max(rms(E(away)), eps);
                end
            end
            info.peak_floor = ratio;
            info.median     = median(ratio(:), 'omitnan');
            info.cut        = obj.param.corr_floor_medfrac * info.median;
            mask = ratio <= info.cut | isnan(ratio);
        end

        function [mask, info] = gate_neighbour(obj, ~, u, v)
        % OUT  mask  ny x nx bool  true = vector disagrees with its neighbourhood
        %      info  struct        on(bool) pivlab_neigh_thresh(float)
        %
        % why  Westerweel-Scarano universal median : a vector against the residuals
        %      of its own neighbourhood, so the threshold means the same thing in a
        %      fast region and a slow one
        % why  reads the FIELD, not the correlation : hence stage 2, after the
        %      three that clean the neighbourhood it is about to read
            info = struct('on', ~isempty(obj.param.pivlab_neigh_thresh), ...
                'pivlab_neigh_thresh', obj.param.pivlab_neigh_thresh);
            mask = false(size(u));
            if ~info.on; return; end

            [filt_u, filt_v] = postproc.PIVlab_postproc(u, v, 1, 1, [], false, ...
                [], true, obj.param.pivlab_neigh_thresh);
            mask = (isnan(filt_u) | isnan(filt_v)) & ~(isnan(u) | isnan(v));
        end
    end
end
