classdef piv_engine < handle
% IN   name-value pairs, any of the settings below
% OUT  obj  configured correlation engine; feed it a stack, get a displacement
%           field. It knows nothing about events, sleep state or the session
%
%   Every PIV setting lives here, so the pipeline stops carrying duplicates of
%   window_sizes and the gate thresholds. Results of the last run are read-only.

    properties  % settings
        % [window step] per pass, coarse to fine. ONE ROW IS ONE PASS: a single
        % [40 20] is a lone 40 px window with no refinement, which leaves the
        % correlation peak too wide to place the displacement
        windows (:,:) double {mustBePositive} = [40 20; 20 10; 12 6]
        exclmask logical = []          % true = EXCLUDE that pixel (PIVlab convention)
        roirect  double  = []          % [x y w h], [] = whole frame
        % Reject a vector shorter than min_snr x its correlation peak half-width.
        % corr2 cannot see this: blur keeps it high while the peak spreads
        min_snr  double = 1            % [] = off
        corr_thr double = 0.5          % corr2 reject threshold
        use_gpu      (1,1) logical = true
        repeat       (1,1) double  = 1
        do_pad       (1,1) double  = 1
        subpixfinder (1,1) double  = 1
        mask_auto    (1,1) double  = 0
        imdeform     (1,:) char    = '*spline'
    end

    properties (SetAccess = private)  % filled by the last run
        uv          % H x W x 2, (:,:,1) = u horizontal, (:,:,2) = v vertical
        corr        % H x W correlation strength, NaN off the vector grid
        planes      % per-pair correlation planes, [] unless kept
        xtable      % vector-grid x coords
        ytable      % vector-grid y coords
        typevector  % 1 = valid, 0 = fully masked
        source      % what the last run was given: stack size, frame indices
    end

    methods
        function obj = piv_engine(opt)
            arguments
                opt.windows      (:,:) double {mustBePositive} = [40 20; 20 10; 12 6]
                opt.exclmask     logical = []
                opt.roirect      double  = []
                opt.min_snr      double  = 1
                opt.corr_thr     double  = 0.5
                opt.use_gpu      (1,1) logical = true
                opt.repeat       (1,1) double  = 1
                opt.do_pad       (1,1) double  = 1
                opt.subpixfinder (1,1) double  = 1
                opt.mask_auto    (1,1) double  = 0
                opt.imdeform     (1,:) char    = '*spline'
            end

            % 1. Checking the pass ladder before anything is stored
            if size(opt.windows, 2) > 2
                error('piv_engine:windowShape', ...
                    'windows must be Px1 or Px2 [window step], got %d columns.', ...
                    size(opt.windows, 2));
            end

            % 2. Copying every option onto the object
            for f = string(fieldnames(opt))'
                obj.(f) = opt.(f);
            end
        end
    end
end
