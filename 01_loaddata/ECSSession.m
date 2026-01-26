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
            obj.img_param.imgstarttime = obj.img_param.imgstartframe / obj.img_param.record_fps;
            obj.img_param.imgendframe = str2double(obj.info.loadend);
            obj.img_param.imgendtime = obj.img_param.imgendframe / obj.img_param.record_fps;
            obj.img_param.pixel2um = str2double(obj.info.objpix);

            % Get group projection and duration
            groupz = str2double(obj.info.groupz);
            fduration_ms = str2double(obj.info.fduration); % Usually in ms

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
            obj.img_param.imgstarttime = (start_frame_idx * fduration_ms) / 1000;

            % End time (seconds)
            % Use valid_raw_frames to exclude deprecated ones
            duration_sec = (valid_raw_frames * fduration_ms) / 1000;
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
                'line_fwhm.mat', 'pax_fwhm'; % Map internal name to class property
                'roilist.mat', 'roilist';
                'polarcluster.mat', 'polarcluster';
                'radon_analysis.mat', 'radon_analysis'; % Use new file name
                'radon_result.mat', 'radon_analysis'; % Fallback for legacy files
                };

            for i = 1:size(files_to_load, 1)
                fname = files_to_load{i, 1};
                prop_name = files_to_load{i, 2};
                fpath = fullfile(savepath, fname);

                if isfile(fpath)
                    tmp = load(fpath);
                    % Try to match property name first
                    if isfield(tmp, prop_name)
                        obj.(prop_name) = tmp.(prop_name);
                    else
                        % Fallback: take the first field
                        fields = fieldnames(tmp);
                        obj.(prop_name) = tmp.(fields{1});
                    end
                end
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
    end

    methods (Access = private)
        function obj = setup_directories(obj, base_path)
            % Integated from manage_directories.m
            % Now merges into dir_struct (inherited from mdfExtractLoader)

            % obj.dir_struct is already initialized by superclass constructor
            obj.dir_struct.load_dir = base_path;
            obj.dir_struct.primary_analysis = fullfile(base_path, 'primary_analysis');

            % Ensure primary_analysis directory exists
            if ~exist(obj.dir_struct.primary_analysis, 'dir')
                mkdir(obj.dir_struct.primary_analysis);
                disp(['Created directory: ', obj.dir_struct.primary_analysis]);
            end

            % Create unique figure directory to preserve results
            timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
            obj.dir_struct.figures = fullfile(obj.dir_struct.primary_analysis, ['figures_', timestamp]);

            if ~exist(obj.dir_struct.figures, 'dir')
                mkdir(obj.dir_struct.figures);
                disp(['Created figure directory: ', obj.dir_struct.figures]);
            end

            % Create subdirectories
            subdirs = {'fwhm', 'radon_figures', 'polarcluster', 'roi_fig'};
            fields  = {'figures_fwhm', 'figures_radon', 'figures_polarcluster', 'figures_roi'};

            for i = 1:numel(subdirs)
                dpath = fullfile(obj.dir_struct.figures, subdirs{i});
                if ~exist(dpath, 'dir')
                    mkdir(dpath);
                end
                obj.dir_struct.(fields{i}) = dpath;
            end
        end
    end
end
