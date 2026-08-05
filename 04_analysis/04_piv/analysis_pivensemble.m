classdef analysis_pivensemble < handle
% IN   imgstack   H x W x T whole recording; only the frames needed are kept
%      frames     [from to], the event endpoints
%      fps        frames per second, for px/frame -> px/s
%      pixel2um   microns per pixel, for px -> um
%      halfwin    frames either side of each endpoint, default 2
%      name-value any pre-processing setting, applied once here
% OUT  obj        configured ensemble-correlation engine. It knows nothing about
%                 events, sleep state or the session -- only two frame numbers
%
%   The margin is the object's business, not the caller's. Correlating from->to
%   needs frames outside [from to], so it takes the whole recording and cuts
%   from-halfwin : to+halfwin itself. The raw stack is never held: what survives
%   is the filtered frames, and only the ones the two schemes actually correlate.
%
%   run() fills BOTH results, because neither answers on its own:
%     between      leading 2*halfwin+1 frames against the trailing ones, so every
%                  pair already spans the whole event
%     consecutive  pairs (1,2) (3,4) ..., each one frame interval, scaled up by
%                  the interval count. Reaching the same place between reports is
%                  what says the move accumulated instead of arriving as a jump
%   BOTH uv fields are the total displacement across the event, in px. The scaling
%   happens in here and the factor applied is recorded as R.scale, so nothing
%   outside ever multiplies a field and no rate can be mistaken for a total.
%
%   between.stack is the interleaved (from_1, to_1, from_2, to_2, ...) frames.
%   Stepped one slice at a time it is a blink comparator, which is the cheapest
%   way to see whether a wall really moved, so it is kept rather than rebuilt.

    properties
        % Correlation settings, flat. A pivlab_ prefix means the value goes
        % straight to PIVlab and is theirs to define; anything unprefixed is
        % defined here, so a call site shows its own provenance. Everything in
        % here is read at run time and can be changed between runs.
        %
        % min_snr is OFF. It rejects a vector shorter than that many correlation
        % peak half-widths, and a short vector is exactly what a piece of tissue
        % that did not move looks like. Stationary structure is a measurement,
        % not a failure, and where the deformation stops is part of the result --
        % so the gate was deleting the answer. Three other measurements point the
        % same way: it is the only filter with a direction preference (with the
        % translation left in, 11 o'clock survived 8% against 4 o'clock's 54%),
        % it removed all 2002 consecutive vectors while corr2 and the local median
        % passed every one, and in an ENSEMBLE the peak width carries the spread
        % of displacements in the window rather than a localisation limit.
        % What remains is corr2 (out-of-plane decorrelation), the
        % Westerweel-Scarano local median (spatial outliers) and track_pct.
        % Set a number to switch min_snr back on for a comparison; between only.
        %
        % track_pct drops the windows with the weakest Shi-Tomasi lambda2, read
        % from the reference image before any correlation. It is the gate min_snr
        % was reaching for and could not express: a window with no corner to lock
        % onto still correlates well, and its vector is placed arbitrarily along
        % the degenerate direction, while a window that simply did not move keeps
        % a strong lambda2 and survives. A percentile rather than an absolute
        % level because lambda2's scale depends on the pre-processing, so no
        % calibrated number carries between sessions. That makes the amount
        % dropped a choice, not a threshold the data sets -- it always removes
        % this fraction, even when every window is fine
        param = struct( ...
            'pivlab_windows',      [40 20; 20 10; 12 6], ... % [window step] PER PASS
            'pivlab_mask_auto',    0, ...
            'pivlab_subpixfinder', 1, ...          % 1 = 3-point gauss, 2 = 2-D fit
            'pivlab_imdeform',     '*spline', ...  % deformation kernel, pass 2 on
            'pivlab_do_pad',       1, ...          % zero-pad -> linear, not circular
            'pivlab_repeat',       1, ...          % pass 1 re-run on 4 quarter shifts
            'pivlab_corr_thr',     0.5, ...        % postproc.PIVlab_correlation_filter
            'pivlab_mask',         [], ...         % true = pixel excluded, [] = none
            'translation_tile',    40, ...  % px tile for the field-wide translation,
                                     ...    % out before the gates. [] leaves it in,
                                     ...    % which lets it pick the surviving sectors
            'min_snr',             [], ...  % OFF; see above. A number switches it on
            'track_pct',           20)      % drop this bottom % of Shi-Tomasi
                                            % lambda2; [] = off. See above
    end

    properties (SetAccess = private)
        % What went in. Settled at construction because each one decides which
        % frames get prepared; changing them afterwards could not rebuild the
        % stacks, so they are not left where they look adjustable
        input = struct( ...
            'fps',      [], ...   % frames per second
            'pixel2um', [], ...   % microns per pixel
            'frames',   [], ...   % [from to] in the recording
            'halfwin',  [], ...   % frames either side of each endpoint
            'span',     [], ...   % the indices actually cut out
            'preproc',  [])       % the pre-processing that was applied

        between     % stack + total displacement across the event; see the header
        consecutive % stack + displacement per frame interval; see the header
    end

    methods
        function obj = analysis_pivensemble(imgstack, frames, fps, pixel2um, opt)
        %   Pre-processing is contrast work only, nothing that moves a pixel:
        %   CLAHE lifts the sparse texture the correlation keys on and Wiener
        %   takes the shot noise with it. It runs here rather than as a method
        %   because it has to happen exactly once, before anything else.
            arguments
                imgstack {mustBeNumeric, mustBeNonempty}
                frames   (1,2) double {mustBePositive, mustBeInteger}
                fps      (1,1) double {mustBePositive}
                pixel2um (1,1) double {mustBePositive}
                opt.halfwin      (1,1) double {mustBeNonnegative} = 2
                opt.clahe        (1,1) logical = true   % lifts the sparse texture
                opt.clahe_size   (1,1) double  = 64     % tile size, px
                opt.wiener       (1,1) logical = true   % denoise; also the low-pass
                opt.wiener_size  (1,1) double  = 3      % px; sigma is half of it
                opt.prctile_low  (1,1) double  = 0.05   % imadjust, on the 0~1 frame
                opt.prctile_high (1,1) double  = 0.95
            end
            % 0. The frames the two schemes need, margin included
            hw = opt.halfwin;
            n  = 2*hw + 1;
            span = frames(1)-hw : frames(2)+hw;
            if span(1) < 1 || span(end) > size(imgstack, 3)
                error('analysis_pivensemble:margin', ...
                    ['event %d-%d needs frames %d-%d for halfwin %d, and the ' ...
                     'recording is %d long. Clipping would move an endpoint, so ' ...
                     'the event is refused instead.'], frames(1), frames(2), ...
                    span(1), span(end), hw, size(imgstack, 3));
            end
            if numel(span) < 2*n
                error('analysis_pivensemble:framesetOverlap', ...
                    ['event %d-%d gives %d frames but halfwin %d needs %d; below ' ...
                     'that the two framesets overlap and a frame correlates with ' ...
                     'itself.'], frames(1), frames(2), numel(span), hw, 2*n);
            end

            obj.input.fps      = fps;
            obj.input.pixel2um = pixel2um;
            obj.input.frames   = frames;
            obj.input.halfwin  = hw;
            obj.input.span     = span;
            obj.input.preproc  = rmfield(opt, 'halfwin');

            % 1. Filter once. The raw stack is not kept, so the two cannot drift
            S = opticalflow_preprocess(imgstack(:, :, span), ...
                clahe       = opt.clahe,       clahe_size   = opt.clahe_size, ...
                wiener      = opt.wiener,      wiener_size  = opt.wiener_size, ...
                prctile_low = opt.prctile_low, prctile_high = opt.prctile_high);

            % 2. Consecutive keeps the span as it stands; pivensemble pairs 1:2:N
            N  = numel(span);
            np = floor(N/2);
            obj.consecutive = struct('stack', S, ...
                'pair_frames', span([(1:2:2*np)', (2:2:2*np)']));

            % 3. Between interleaves the two ends into (a_1, b_1, a_2, b_2, ...)
            a = 1:n;
            b = N-n+1 : N;
            [H, W] = size(S, 1, 2);
            inter = zeros(H, W, 2*n);
            inter(:, :, 1:2:end) = S(:, :, a);
            inter(:, :, 2:2:end) = S(:, :, b);
            obj.between = struct('stack', inter, 'pair_frames', span([a(:), b(:)]));
        end

        function run(obj)
        % OUT  none; fills both results with uv, uv_ungated, corr, planes,
        %      common_mode, leaving stack and pair_frames as they were
        %
        %   Both, always. between measures the displacement and consecutive says
        %   whether it accumulated: a jump would put all but one pair plane at
        %   zero, and summing planes is a majority vote, so the jump would be
        %   voted down to nothing. A consecutive sum that reaches between's answer
        %   is the evidence that it did not teleport.
        %
        %   The width gate goes to between only. Consecutive measures a per-interval
        %   rate that is 0.02 px by construction, so a threshold in units of peak
        %   width has nothing to divide there -- see param.
            obj.between     = obj.correlate(obj.between,     obj.param.min_snr, 1);
            % Scaled by the intervals between the endpoints, not by the pair count:
            % the pairs sample every other interval, so counting pairs halves it
            obj.consecutive = obj.correlate(obj.consecutive, [], ...
                diff(obj.input.frames));
        end
    end

    methods (Static)
        function [E, info] = plane_at(R, x, y)
        % IN   R      a result struct, obj.between or obj.consecutive
        %      x, y   image coordinates in px; the nearest grid point is used
        % OUT  E      IA x IA ensembled plane, centre cell = zero displacement
        %      info   grid indices, the coordinates used, uv and validity
        %
        %   Takes the result rather than picking one, so the caller names which
        %   field it is reading and the class carries no selector.
        %
        %   The planes were captured before validation, so a window the gates
        %   threw away still has one. That is the point of looking at it.
            if ~isfield(R, 'planes') || isempty(R.planes)
                error('analysis_pivensemble:noPlanes', 'run() has not run.');
            end
            cp = R.planes;
            d  = hypot(cp.xtable - x, cp.ytable - y);
            [~, k] = min(d(:));
            [r, c] = ind2sub(size(d), k);
            E = double(sum(cp.maps(:, :, r, c, :), 5));
            if nargout > 1
                ri = round(cp.ytable(r, c));   ci = round(cp.xtable(r, c));
                info = struct('row', r, 'col', c, ...
                    'x', cp.xtable(r, c), 'y', cp.ytable(r, c), ...
                    'uv', squeeze(R.uv(ri, ci, :))', ...
                    'valid', cp.typevector(r, c) == 1);
            end
        end
    end

    methods (Access = private)
        function R = correlate(obj, R, min_snr, scale)
        % IN   R        a result struct carrying stack and pair_frames
        %      min_snr  width-gate threshold for this result, [] = skip that gate
        %      scale    factor onto the displacement, so every result comes out as
        %               the total across the event whatever the pairing was
        % OUT  R        the same, with uv, uv_ungated, corr, planes, common_mode,
        %               scale
        %
        %   The order is the whole reason this lives here instead of inside
        %   corr_ensemble: the translation comes out BEFORE a gate that keys off
        %   vector length, or that gate keeps whichever sectors point along it.
        %   Measured on ev2 with the translation left in, 11 o'clock survived 8%
        %   against 4 o'clock's 54%.
            p = obj.param;
            [H, W, ~] = size(R.stack);

            % 1. Correlate only; corr_ensemble's own gates stay off
            [~, ~, cp] = corr_ensemble(R.stack, ...
                window_sizes = p.pivlab_windows,      exclmask  = p.pivlab_mask, ...
                subpixfinder = p.pivlab_subpixfinder, mask_auto = p.pivlab_mask_auto, ...
                imdeform     = p.pivlab_imdeform,     repeat    = p.pivlab_repeat, ...
                do_pad       = p.pivlab_do_pad,       validate  = false, ...
                save_corrmaps = true);

            % 2. Field-wide translation out
            [u0, v0, cm] = obj.strip_translation(cp.utable, cp.vtable, ...
                cp.xtable(1,2) - cp.xtable(1,1));

            % 3. Untrackable windows out first, read from the images so the
            %    judgement owes nothing to the correlation it is judging. BOTH
            %    ends: a correspondence needs texture at the start and at the
            %    finish, so a window that had a corner and lost it by the second
            %    frame is no more trackable than one that never had one. Gate on
            %    the weaker of the two; lambda2 keeps them separately, (:,:,1) =
            %    leading frames, (:,:,2) = trailing, so which end failed is visible
            u = u0;   v = v0;
            R.lambda2 = [];
            if ~isempty(p.track_pct)
                win = p.pivlab_windows(end, 1);
                R.lambda2 = cat(3, ...
                    piv_trackability(mean(R.stack(:,:,1:2:end), 3), cp.xtable, ...
                        cp.ytable, win, p.pivlab_mask), ...
                    piv_trackability(mean(R.stack(:,:,2:2:end), 3), cp.xtable, ...
                        cp.ytable, win, p.pivlab_mask));
                worst = min(R.lambda2, [], 3);
                weak  = worst <= prctile(worst(:), p.track_pct);
                u(weak) = NaN;   v(weak) = NaN;
            end

            % 4. corr2 gate, then the width gate if this result gets one, then the
            %    local median. The gates run on the unscaled field, because the
            %    width threshold is in px of what was actually correlated
            pw = [];
            if ~isempty(min_snr); pw = piv_peak_width(cp.maps); end
            [u, v] = piv_validate(u, v, cp.corr, corr_thr = p.pivlab_corr_thr, ...
                peak_width = pw, min_snr = min_snr, neigh_thresh = 2);

            % 5. Stamp both fields, scaled to the event. uv_ungated is the field
            %    as correlated, so the difference is exactly what the gates took
            R.uv_ungated  = obj.stamp(u0, v0, cp, H, W) * scale;
            R.uv          = obj.stamp(u,  v,  cp, H, W) * scale;
            R.common_mode = cm * scale;
            R.scale       = scale;
            R.planes      = cp;                 % planes stay unscaled, as correlated

            % 6. Correlation strength where a vector survived
            R.corr = NaN(H, W);
            keep = (cp.typevector == 1) & ~isnan(u) & ~isnan(v);
            rows = max(1, min(H, round(cp.ytable(keep))));
            cols = max(1, min(W, round(cp.xtable(keep))));
            R.corr(sub2ind([H, W], rows, cols)) = cp.corr(keep);
        end

        function M = stamp(~, u, v, cp, H, W)
        % OUT  M  H x W x 2 with each vector at its own pixel, NaN elsewhere
        %
        %   No interpolation. A vector belongs to the pixel its window centre sits
        %   on, and the gaps between stay NaN so nothing invented goes downstream.
            M = NaN(H, W, 2);
            keep = (cp.typevector == 1) & ~isnan(u) & ~isnan(v);
            rows = max(1, min(H, round(cp.ytable(keep))));
            cols = max(1, min(W, round(cp.xtable(keep))));
            idx  = sub2ind([H, W], rows, cols);
            M(idx)       = u(keep);
            M(idx + H*W) = v(keep);
        end

        function [u, v, cm] = strip_translation(obj, u, v, step)
        % IN   u, v  raw sparse field on the vector grid
        %      step  grid spacing in px, so translation_tile can be read in cells
        % OUT  u, v  the same field with the translation taken out
        %      cm    [u v] that was removed; [0 0] when translation_tile is []
        %
        %   A radially symmetric dilation cancels to zero, so what is left over the
        %   whole field is a rigid shift: registration drift, or the tissue really
        %   translating. Those two are not separable from one event, which is why
        %   cm is kept rather than discarded.
        %
        %   Tiled, not a plain median. Vector density runs 2-7x uneven across the
        %   frame and tracks local displacement, so a plain median leans toward
        %   wherever the field is strongest. One median per tile, then the median
        %   of tiles, weights by area instead.
            cm = [0 0];
            if isempty(obj.param.translation_tile); return; end

            ok = ~isnan(u) & ~isnan(v);
            % 1. Tile index of every grid point
            tile = max(1, round(obj.param.translation_tile / step));
            [ny, nx] = size(u);
            ti = ceil((1:ny)'/tile) * ones(1, nx);
            tj = ones(ny, 1) * ceil((1:nx)/tile);
            id = (ti - 1)*max(tj(:)) + tj;
            % 2. One median vector per tile; tiles too sparse to have one drop out
            cnt = accumarray(id(ok), 1);
            mu  = accumarray(id(ok), u(ok), [], @median);
            mv  = accumarray(id(ok), v(ok), [], @median);
            use = cnt >= 20;
            if ~any(use); return; end

            cm = [median(mu(use)), median(mv(use))];
            u  = u - cm(1);
            v  = v - cm(2);
        end
    end
end
