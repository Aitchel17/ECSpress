classdef analysis_piv < handle
    %ANALYSIS_PIV  PIV / optical-flow analysis of ONE image stack.
    %   The stack is held by reference (copy-on-write), not duplicated. run_ensemble
    %   and run_wavelet both write displacement_map, stamped with .method, .units and
    %   the parameters that produced it; add_divergence consumes whichever ran last.
    %
    % IN   imgstack        H x W x N stack, already preprocessed
    %      frame_interval  seconds per frame (dt); scales velocity and t_axis
    %      resolution      microns per pixel; scales displacement
    % OUT  obj             analysis_piv handle

    properties
        imgstack                 % H x W x N working stack (reference, no copy)
        frame_interval           % seconds per frame (dt)
        resolution               % microns per pixel

        vectormap = struct( ...
            'displacement_map', [], ...  % dense field, um/s: .data .method .units .param
            'divergent_map',    [], ...  % dense divergence, 1/s: .data .units .param
            'corr_map',         [], ...  % ensemble confidence: .data .units .param
            'tracks',           [])      % sparse KLT tracks, px & frame: .data .method .units .param

        % per-algorithm options; set before calling the matching run_*
        param = struct( ...
            'ensemble',   struct('window_sizes', [64 32; 32 16; 16 8], 'subpixfinder', 1, ...
                                 'mask_auto', 0, 'imdeform', '*spline', 'repeat', 0, ...
                                 'do_pad', 1, 'use_gpu', false), ...
            'wavelet',    struct('pyramid_levels', 3, 'smoothness', 50), ...
            'klt',        struct('targetN', 50, 'patchR', 15, 'nccThr', 0.5, 'scoreThr', 0.4, ...
                                 'maxDisp', 12, 'minReseedDist', 15, 'minEigenQuality', 0.05, ...
                                 'maxBidiErr', 2, 'blockSize', [31 31], 'pyrLevels', 3), ...
            'validation', struct('corr_thr', 0.5, 'neigh_thresh', 2), ...
            'roi',        struct('roirect', []), ...
            'divergence', struct('mask', []))
    end

    methods
        function obj = analysis_piv(imgstack, frame_interval, resolution)
            %ANALYSIS_PIV  Store the stack and its time and space scales.
            arguments
                imgstack       {mustBeNumeric, mustBeNonempty}
                frame_interval (1,1) double
                resolution     (1,1) double
            end
            obj.imgstack       = imgstack;
            obj.frame_interval = frame_interval;
            obj.resolution     = resolution;
        end

        function run_ensemble(obj)
            %RUN_ENSEMBLE  Ensemble multipass cross-correlation PIV.
            %   Writes displacement_map (dense, time-averaged, um/s) and corr_map.
            p = obj.param.ensemble;
            r = obj.param.roi;
            % piv_corr_ensemble correlates and nothing else, so the gating is
            % spelled out here where param.validation already claimed to describe
            % it. Order is PIVlab's: correlation filter, then the local median
            cp = piv_corr_ensemble(obj.imgstack, ...
                window_sizes = p.window_sizes, roirect = r.roirect, ...
                subpixfinder = p.subpixfinder, mask_auto = p.mask_auto, ...
                imdeform = p.imdeform, repeat = p.repeat, do_pad = p.do_pad, ...
                use_gpu = p.use_gpu);
            q = obj.param.validation;
            [cp.utable, cp.vtable] = piv_validate(cp.utable, cp.vtable, cp.corr, ...
                corr_thr = q.corr_thr, neigh_thresh = q.neigh_thresh);
            [uv, cm] = vfield_stamp(cp);
            uv = uv * (obj.resolution / obj.frame_interval);   % px/frame -> um/s
            obj.vectormap.displacement_map = struct('data', uv, 'method', 'ensemble', 'units', 'um/s', 'param', p);
            obj.vectormap.corr_map         = struct('data', cm, 'units', 'corr', 'param', obj.param.validation);
        end

        function run_wavelet(obj)
            %RUN_WAVELET  Wavelet-based optical flow.
            %   Writes displacement_map (dense per-frame, H x W x T x 2, um/s).
            p = obj.param.wavelet;
            r = obj.param.roi;
            flow = opticalflow_wlet(obj.imgstack, ...
                pyramid_levels = p.pyramid_levels, smoothness = p.smoothness, roirect = r.roirect);
            flow = flow * (obj.resolution / obj.frame_interval);   % px/frame -> um/s
            obj.vectormap.displacement_map = struct('data', flow, 'method', 'wavelet', 'units', 'um/s', 'param', p);
        end

        function run_klt(obj)
            %RUN_KLT  KLT feature tracking with re-seeding.
            %   Writes tracks (sparse points); left in pixels and frames, unscaled.
            p = obj.param.klt;
            [tracks, ~] = klt_track_reseed(obj.imgstack, ...
                targetN = p.targetN, patchR = p.patchR, nccThr = p.nccThr, ...
                scoreThr = p.scoreThr, maxDisp = p.maxDisp, minReseedDist = p.minReseedDist, ...
                minEigenQuality = p.minEigenQuality, maxBidiErr = p.maxBidiErr, ...
                blockSize = p.blockSize, pyrLevels = p.pyrLevels);
            obj.vectormap.tracks = struct('data', tracks, 'method', 'klt', 'units', 'px,frame', 'param', p);
        end

        function add_divergence(obj, mask)
            %ADD_DIVERGENCE  Divergence of the dense displacement field, in 1/s.
            %   Run run_ensemble or run_wavelet first; mask defaults to param.divergence.mask.
            arguments
                obj
                mask = obj.param.divergence.mask
            end
            if isempty(obj.vectormap.displacement_map)
                error('analysis_piv:noField', 'Run ensemble() or wavelet() before add_divergence().');
            end
            % field is um/s on a pixel grid, so divide by resolution (um/px) -> 1/s
            div = vfield_divergence_fd(obj.vectormap.displacement_map.data, mask = mask) / obj.resolution;
            obj.vectormap.divergent_map = struct('data', div, 'units', '1/s', 'param', struct('mask', mask));
            obj.param.divergence.mask = mask;
        end

        function t = t_axis(obj)
            %T_AXIS  Frame times in seconds, from frame_interval.
            t = (0:size(obj.imgstack, 3) - 1) * obj.frame_interval;
        end

        function save2disk(obj, name, savepath)
            %SAVE2DISK  Save the object to savepath/name.mat as variable analysis_piv.
            analysis_piv = obj;
            save(fullfile(savepath, [name, '.mat']), 'analysis_piv');
        end
    end
end
