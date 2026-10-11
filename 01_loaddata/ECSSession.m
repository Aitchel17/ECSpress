classdef ECSSession < mdfExtractLoader
    % ECSSession Specialized loader for ECSPRESS image processing pipeline
    %   Inherits from mdfExtractLoader to provide specific context for
    %   Two-Photon image processing, including parameter calculation and
    %   primary result loading.

    properties
        img_param       % Structure containing imaging parameters (fps, um/pixel, time axis)

        % directories     % Removed: merged into dir_struct (inherited from mdfExtractLoader)
        stackch1        % Image stack channel 1
        stackch2        % Image stack channel 2

        % Primary Analysis Objects
        roilist
        pax_fwhm
        polarcluster
        radon_analysis

        % Peripheral products and tissue motion, by load_peripheral and load_motion; [] when absent
        primary_analog  % analysis_analog | []   peripheral/analysis_analog.mat
        camera_session  % analysis_camera | []   peripheral/analysis_camera.mat
        sleep_score     % 1x1 struct | []        peripheral/sleep_score.mat, as load returns it
        motion          % 1x1 struct | []        t_s, rostral_um, lateral_um : N x 1 float each
    end

    methods
        function obj = ECSSession(extract_dir)
            % ECSSession Construct an instance of this class
            %   extract_dir: Path to the mdfExtracted folder

            % Superclass constructor must be called unconditionally
            if nargin == 0
                args = {};
            else
                args = {extract_dir};
            end
            obj@mdfExtractLoader(args{:});

            % Setup Directories
            if nargin > 0
                obj = obj.setup_directories(extract_dir);
            end

            % Calculate Imaging Parameters
            obj = obj.calculate_img_params();
        end

        function obj = calculate_img_params(obj)
            % Logic ported from load_mdfextract.m
            obj.img_param = struct();
            obj.img_param.save_fps = str2double(obj.info.savefps);
            obj.img_param.record_fps = str2double(obj.info.fps);
            obj.img_param.imgstartframe = str2double(obj.info.loadstart);
            obj.img_param.imgendframe = str2double(obj.info.loadend);
            objpix = strsplit(obj.info.objpix);
            obj.img_param.pixel2um = str2double(objpix{1});
            

            % Get group projection and duration
            groupz = str2double(obj.info.groupz);
            fduration_sec = str2double(obj.info.fduration(1:end-1)); % Usually in ms

            % Frames
            start_frame_idx = str2double(obj.info.loadstart);
            end_frame_idx = str2double(obj.info.loadend);
            total_raw_frames = end_frame_idx - start_frame_idx + 1;

            % Deprecated frames at the end
            deprecated_frames = mod(total_raw_frames, groupz);
            valid_raw_frames = total_raw_frames - deprecated_frames;

            % Resulting frames in stack (after grouping)
            num_stack_frames = valid_raw_frames / groupz;

            % Time calculation (fduration is usually ms per frame)
            % Start time (seconds)
            obj.img_param.imgstarttime = (start_frame_idx * fduration_sec);

            % End time (seconds)
            % Use valid_raw_frames to exclude deprecated ones
            duration_sec = (valid_raw_frames * fduration_sec);
            obj.img_param.imgendtime = obj.img_param.imgstarttime + duration_sec;

            % T-axis
            % Time points corresponding to the CENTER of each group? Or Start?
            % Usually linspace from start to end.
            obj.img_param.taxis = linspace(obj.img_param.imgstarttime, obj.img_param.imgendtime, num_stack_frames);
        end


        function obj = load_primary_results(obj, extractfolder_path)
            % Logic ported from load_primaryresult.m
            % Fixes bug: uses 'obj' (self) instead of undefined 'mdfextract'

            if nargin < 2
                % Default to the folder of the MDF data if not provided
                extractfolder_path = obj.dir_struct.info; % This is a file path, we need folder
                [extractfolder_path, ~, ~] = fileparts(extractfolder_path);
            end

            savepath = fullfile(extractfolder_path, 'primary_analysis');
            if ~exist(savepath, 'dir')
                mkdir(savepath);
            end

            % --- Load Other Results ---
            files_to_load = {
                'paxfwhm.mat', 'line_fwhm'; % Map internal name to class property
                'roilist.mat', 'roilist';
                'polarcluster.mat', 'polarcluster';
                };

            for i = 1:size(files_to_load, 1)
                fname = files_to_load{i, 1};
                prop_name = files_to_load{i, 2};
                fpath = fullfile(savepath, fname);

                if isfile(fpath)
                    tmp = load(fpath);
                    % Try to match property name first
                    if isfield(tmp, prop_name)
                        if strcmp(fname, 'paxfwhm.mat')
                            obj.pax_fwhm = tmp.(prop_name);
                        else
                            obj.(prop_name) = tmp.(prop_name);
                        end

                    else
                        % Fallback: take the first field
                        fields = fieldnames(tmp);
                        obj.(prop_name) = tmp.(fields{1});
                    end
                end
            end
            % A polarcluster.mat written before analysis_clusterpolar was a class holds a struct
            if isstruct(obj.polarcluster)
                obj.polarcluster = analysis_clusterpolar.from_struct(obj.polarcluster);
            end

            % Radon
            % Try loading the new object-based file first
            if isfile(fullfile(savepath, 'radon_analysis.mat'))
                tmp = load(fullfile(savepath, 'radon_analysis.mat'));
                obj.radon_analysis = tmp.radon_analysis;
                % Fallback to the old struct-based file
            elseif isfile(fullfile(savepath, 'radon_result.mat'))
                tmp = load(fullfile(savepath, 'radon_result.mat'));
                obj.radon_analysis = tmp.radon_result;
            end

            % Ensure roilist is initialized so downstream code (addormodifyroi) doesn't crash
            if isempty(obj.roilist)
                try
                    % roi_handle expects the full path to the .mat file
                    obj.roilist = roi_handle(fullfile(savepath, 'roilist.mat'));
                catch
                    warning('Could not initialize default roi_handle.');
                end
            end
        end

        function obj = load_peripheral(obj)
            % The three products peripheral_batch and the sleep scoring write; one left [] when its file is absent
            peripheral_dir = fullfile(obj.dir_struct.load_dir, 'peripheral');
            product = { ...              % file, variable inside it ('' = the whole file), property
                'analysis_analog.mat', 'primary_analog', 'primary_analog'; ...
                'analysis_camera.mat', 'camera_session', 'camera_session'; ...
                'sleep_score.mat',     '',               'sleep_score'};
            for k = 1:size(product, 1)
                file_path = fullfile(peripheral_dir, product{k, 1});
                if ~isfile(file_path)
                    fprintf('   no %s\n', file_path);
                    continue
                end
                if isempty(product{k, 2})
                    obj.(product{k, 3}) = load(file_path);
                else
                    obj.(product{k, 3}) = load(file_path, product{k, 2}).(product{k, 2});
                end
            end
        end

        function obj = load_motion(obj)
            % Tissue displacement from <session>_motion.txt, um from the session median : rostral = image top (-row), lateral = image right (+col)
            motion_file = dir(fullfile(obj.dir_struct.load_dir, '*_motion.txt'));
            if isempty(motion_file)
                fprintf('   no *_motion.txt in %s\n', obj.dir_struct.load_dir);
                return
            end
            motion_path = fullfile(motion_file(1).folder, motion_file(1).name);
            [disp_row, disp_col, ~, drift_fps] = io_readmotion(motion_path);
            frame_idx = (0:numel(disp_row) - 1)';
            row_shift = disp_row(:) - median(disp_row);
            col_shift = disp_col(:) - median(disp_col);
            obj.motion = struct( ...
                't_s', frame_idx / drift_fps, ...
                'rostral_um', -row_shift * obj.img_param.pixel2um, ...
                'lateral_um', col_shift * obj.img_param.pixel2um);
        end
    end

    methods (Access = private)
        function obj = setup_directories(obj, base_path)
            obj.dir_struct.load_dir = base_path;
            obj.dir_struct.primary_analysis = fullfile(base_path, 'primary_analysis');

            % Ensure primary_analysis directory exists
            if ~exist(obj.dir_struct.primary_analysis, 'dir')
                mkdir(obj.dir_struct.primary_analysis);
                disp(['Created directory: ', obj.dir_struct.primary_analysis]);
            end


            %% figure directory structure
            % Flat, lazily-created figure folders directly under primary_analysis:
            %   <yyMMdd_HHmmss>_<category>     e.g. 260717_143022_polarcluster
            % The folders are NOT created here. make_fig.save2svg creates each one
            % only when a figure is actually saved into it, so categories (and whole
            % session loads) that produce no figures never leave an empty folder.
            % All categories from one run share the timestamp prefix, so they still
            % sort together.
            delete_emptydir(fullfile(obj.dir_struct.primary_analysis, 'figures_*'));  % clean up legacy nested run folders

            timestamp = char(datetime('now', 'Format', 'yyMMdd_HHmmss'));
            obj.dir_struct.figures = obj.dir_struct.primary_analysis;  % base (kept for back-compat)

            subdirs = {'fwhm', 'radon_figures', 'polarcluster', 'roi_fig'};
            fields  = {'figures_fwhm', 'figures_radon', 'figures_polarcluster', 'figures_roi'};
            for i = 1:numel(subdirs)
                obj.dir_struct.(fields{i}) = fullfile(obj.dir_struct.primary_analysis, [timestamp, '_', subdirs{i}]);
            end
        end
    end
end
