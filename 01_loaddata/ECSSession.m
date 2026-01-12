classdef ECSSession < mdfExtractLoader
    % ECSSession Specialized loader for ECSPRESS image processing pipeline
    %   Inherits from mdfExtractLoader to provide specific context for
    %   Two-Photon image processing, including parameter calculation and
    %   primary result loading.

    properties
        img_param       % Structure containing imaging parameters (fps, um/pixel, time axis)
        results         % Structure containing analysis results (analog, roilist, etc.)
        directories     % Structure containing relevant directory paths
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

            % Initialize container structs
            obj.results = struct();

            % Setup Directories
            if nargin > 0
                obj = obj.setup_directories(extract_dir);
            end

            % Calculate Imaging Parameters
            obj = obj.calculate_img_params();

            % Load Stacks
            if nargin > 0
                obj.stackch1 = obj.loadstack('ch1');
                obj.stackch2 = obj.loadstack('ch2');
            end
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
            obj.img_param.pixel2um = str2double(obj.info.objpix(1:end-2));

            % Note: taxis requires stack size. We can infer it from info or load it.
            % load_mdfextract loaded the stack to get size(stackch1,3).
            % To avoid loading stacks eagerly just for taxis, we can use info.loadend - info.loadstart?
            % Or just calculate it when stacks are actually loaded.
            % For now, let's keep it consistent: We might need to load stack info without full data.
            % But mdfExtractLoader doesn't expose frame count easily without loading.
            % Let's use the frame range from info for now.
            num_frames = obj.img_param.imgendframe - obj.img_param.imgstartframe + 1; % Approximation

            % However, load_mdfextract did: linspace(start, end, size(stack,3))
            % Let's rely on info for T-axis generation to be lazy (efficient).
            obj.img_param.taxis = linspace(obj.img_param.imgstarttime, obj.img_param.imgendtime, num_frames);
        end

        function obj = load_primary_results(obj, extractfolder_path)
            % Logic ported from load_primaryresult.m
            % Fixes bug: uses 'obj' (self) instead of undefined 'mdfextract'

            if nargin < 2
                % Default to the folder of the MDF data if not provided
                extractfolder_path = obj.path_struct.info; % This is a file path, we need folder
                [extractfolder_path, ~, ~] = fileparts(extractfolder_path);
            end

            savepath = fullfile(extractfolder_path, 'primary_analysis');
            if ~exist(savepath, 'dir')
                mkdir(savepath);
            end

            % --- Load Analog ---
            analog_file = fullfile(savepath, 'analog.mat');
            if isfile(analog_file)
                tmp = load(analog_file);
                obj.results.analog = tmp.primary_analog;
            else
                fprintf('analog file does not exist, generating new one...\n');
                % Use INHERITED method
                raw_analog = obj.loadanalog();

                % Create primary_analog object (assuming analysis_analog is available)
                primary_analog = analysis_analog(raw_analog.info, raw_analog.data);

                if isfield(primary_analog.data, 'raw_Air_puff1')
                    primary_analog.airtable = primary_analog.get_airtable('raw_Air_puff1');
                end
                if isfield(primary_analog.data, 'raw_ECoG')
                    primary_analog.ecogspectrum = primary_analog.get_ecogspectrum('raw_ECoG');
                end

                % Save it
                save(analog_file, "primary_analog");
                obj.results.analog = primary_analog;
            end

            % --- Load Other Results ---
            files_to_load = {
                'line_fwhm.mat', 'pax_fwhm'; % Map internal name to class property
                'roilist.mat', 'roilist';
                'polarcluster.mat', 'polarcluster';
                'radon_result.mat', 'radon_analysis'; % Assuming filename
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

            % Manual Polar structs (Storing in results struct for now as they are less standard)
            if isfile(fullfile(savepath, 'manual_polar.mat'))
                tmp = load(fullfile(savepath, 'manual_polar.mat'));
                obj.results.manual_polar_roilist = tmp.roilist;
            end
            if isfile(fullfile(savepath, 'manual_polarstruct.mat'))
                tmp = load(fullfile(savepath, 'manual_polarstruct.mat'));
                obj.results.manual_polarstruct = tmp.manual_polarstruct;
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

            obj.directories = struct();
            obj.directories.load_dir = base_path;
            obj.directories.primary_analysis = fullfile(base_path, 'primary_analysis');

            % Ensure primary_analysis directory exists
            if ~exist(obj.directories.primary_analysis, 'dir')
                mkdir(obj.directories.primary_analysis);
                disp(['Created directory: ', obj.directories.primary_analysis]);
            end

            % Create unique figure directory to preserve results
            timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
            obj.directories.figures = fullfile(obj.directories.primary_analysis, ['figures_', timestamp]);

            if ~exist(obj.directories.figures, 'dir')
                mkdir(obj.directories.figures);
                disp(['Created figure directory: ', obj.directories.figures]);
            end

            % Create subdirectories
            subdirs = {'fwhm', 'radon_figures', 'polarcluster', 'roi_fig'};
            fields  = {'figures_fwhm', 'figures_radon', 'figures_polarcluster', 'figures_roi'};

            for i = 1:numel(subdirs)
                dpath = fullfile(obj.directories.figures, subdirs{i});
                if ~exist(dpath, 'dir')
                    mkdir(dpath);
                end
                obj.directories.(fields{i}) = dpath;
            end
        end
    end
end
