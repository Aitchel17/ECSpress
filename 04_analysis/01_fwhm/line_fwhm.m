classdef line_fwhm < handle
    %LINE_FWHM Summary of this class goes here
    %   Detailed explanation goes here

    properties
        rotatecrop = struct()
        kymograph = struct()
        thickness = struct()
        displacement = struct()
        up_thicker      % bool   trapz(uppvs - downpvs) > 0. WHICH SIDE IS THICKER
                        %        on average, nothing about which one moves. It
                        %        was called updynamic; see CLAUDE_LOG.md
        mask
        idx
        t_axis
        param = struct('input_size',[-1,-1,-1],...
            'line_info',[-1,-1;-1,-1;-1,-999],...
            'idxclean_noisethr', 5, ...
            'idxclean_refmedian', 7); % 3x2 double, {x1,y1; x2,y2; linewidth,-999}
    end


    methods
        function obj = line_fwhm(line_info)
            % Constructor function get roi info, ex: (x1,y1;x2,y2;,thickness,-999)
            arguments
                line_info (3,2)
            end
            obj.param.line_info = line_info;
        end

        function addkymograph(obj,stack_name,stack,mode)
            arguments
                obj
                stack_name (1,1) string {mustBeMember(stack_name,{'lumen','wall', 'pvs', 'outside'})}
                stack
                mode (1,1) string {mustBeMember(mode,{'mean', 'median','max', 'min'})}
            end
            obj.param.input_size = size(stack);
            name_rotatecrop = strcat("rc_",stack_name);
            name_kymograph = strcat("kgph_",stack_name);
            name_normkymograph = strcat("normkgph_",stack_name);
            rotatecroped_stack = analyze_affine_rotate(stack,obj.param.line_info(1:2,:), obj.param.line_info(3,1));
            if strcmp(mode,'mean')
                rawkymograph = squeeze(sum(rotatecroped_stack, 1));
            elseif strcmp(mode,'max')
                rawkymograph = squeeze(max(rotatecroped_stack, [], 1));
            elseif strcmp(mode,'min')
                disp('projection mode min')
                rawkymograph = squeeze(min(rotatecroped_stack,[],1));
            elseif strcmp(mode,'median')
                rawkymograph = squeeze(median(rotatecroped_stack, 1));
            end
            obj.param.(strcat('kymograph_mode_',stack_name)) = mode; 
            obj.rotatecrop.(name_rotatecrop) = rotatecroped_stack;
            obj.kymograph.(name_kymograph) = rawkymograph;
            obj.kymograph.(name_normkymograph) = (rawkymograph-min(rawkymograph,[],1))./max(rawkymograph,[],1); % just for reference to see intermediate step, no offset
        end

        function kymograph_afterprocess(obj,stack_name,window_size)
            arguments
                obj
                stack_name (1,1) string {mustBeMember(stack_name,{'lumen','wall', 'pvs', 'outside'})}
                window_size   (1,2) double = [1 3] % default halfmax
            end
            % parameter:
            low_cutoff = 10;
            high_cutoff = 95;
            medfilt_window = window_size;
            %
            name_kymograph = strcat("kgph_",stack_name);
            name_processed_kymograph = strcat(name_kymograph,'_processed');
            processed_kymograph = medfilt2(obj.kymograph.(name_kymograph),window_size);
            prctile_low = prctile(processed_kymograph,low_cutoff,1);
            prctile_high = prctile(processed_kymograph,high_cutoff,1);
            processed_kymograph = min(processed_kymograph, prctile_high);
            processed_kymograph = max(processed_kymograph, prctile_low);
            processed_kymograph = processed_kymograph - min(processed_kymograph,[],1);
            processed_kymograph = processed_kymograph./max(processed_kymograph,[],1);
            obj.kymograph.(name_processed_kymograph) = processed_kymograph;
            obj.param.(strcat('kymograph_cutoff',stack_name)) = [low_cutoff , high_cutoff];
            obj.param.(strcat('kymograph_medfiltwindow',stack_name)) = medfilt_window;
        end


        function fwhm(obj,stack_name)
            arguments
                obj
                stack_name (1,1) string {mustBeMember(stack_name,{'lumen','wall', 'pvs', 'outside'})}
            end
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            name_kymograph = strcat("kgph_",stack_name,"_processed");
            [tmp.idx, tmp.kgph_mask,tmp.param] = get_bvoutter(obj.kymograph.(name_kymograph));
            % obj.param.bv_thr = threshold;
            % obj.param.bv_offset = offset;
            obj.idx = obj.mergestruct(obj.idx, tmp.idx);
            obj.mask = obj.mergestruct(obj.mask, tmp.kgph_mask);
        end

        function pvsanalysis(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            [tmp.idx, tmp.kgph_mask] = get_pvsoutter(obj.kymograph.kgph_pvs_processed, obj.idx.upperBVboundary, obj.idx.lowerBVboundary);
            obj.idx = obj.mergestruct(obj.idx, tmp.idx);
            obj.mask = obj.mergestruct(obj.mask, tmp.kgph_mask);
        end
        function pvsanalysis_inverted(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            [tmp.idx, tmp.kgph_mask] = get_invertedpvsoutter(obj.kymograph.kgph_pvs_processed, obj.idx.upperBVboundary, obj.idx.lowerBVboundary);
            obj.idx = obj.mergestruct(obj.idx, tmp.idx);
            obj.mask = obj.mergestruct(obj.mask, tmp.kgph_mask);
        end
        function clean_outlier(obj,overwrite)
            arguments
                obj
                overwrite = false
            end
            obj.idx = clean_idxoutlier(obj.idx, overwrite,obj.param.idxclean_refmedian,obj.param.idxclean_noisethr);
        end

        function getdiameter(obj)
            obj.thickness = struct();
        % STORED ANATOMICALLY, on purpose. This used to sort the two sides by
        % thickness and store them as dynamic_pvs / static_pvs, so a claim about
        % MOTION was baked into the field names by a criterion that has no time
        % in it -- trapz integrates the time axis away. Measured over 119
        % sessions the label agreed with its own thickness rule 99.2% of the
        % time and with a motion criterion 65.5%. See CLAUDE_LOG.md and
        % FINDINGS.md. up_thicker still records which side is thicker, so any
        % caller that wants the sorted view can build it and has to say so.
            obj.thickness.bv = obj.idx.clean_lowerBVboundary - obj.idx.clean_upperBVboundary;
            obj.thickness.uppvs = obj.idx.clean_upperBVboundary - obj.idx.clean_pvsupedge_idx;
            obj.thickness.downpvs = obj.idx.clean_pvsdownedge_idx - obj.idx.clean_lowerBVboundary;
            obj.thickness.totalpvs = obj.thickness.uppvs + obj.thickness.downpvs;
            obj.thickness.eps = obj.thickness.totalpvs + obj.thickness.bv;

            % which side is thicker, as its own recorded fact
            difference_pvs = obj.thickness.uppvs - obj.thickness.downpvs;
            difference_pvs = medfilt1(difference_pvs,11); % smoothing the pvs thickness difference
            % area under the curve, so one frame cannot decide it
            obj.up_thicker = trapz(difference_pvs) > 0;

            %%
            obj.thickness.median_totalpvs = median(obj.thickness.totalpvs);
            obj.thickness.median_uppvs = median(obj.thickness.uppvs);
            obj.thickness.median_downpvs = median(obj.thickness.downpvs);
            %%
            obj.thickness.median_bv = median(obj.thickness.bv);
            obj.thickness.bvchanges = obj.thickness.bv - obj.thickness.median_bv;
            obj.thickness.pvschanges_total = obj.thickness.totalpvs - obj.thickness.median_totalpvs;
            obj.thickness.pvschanges_up = obj.thickness.uppvs - obj.thickness.median_uppvs;
            obj.thickness.pvschanges_down = obj.thickness.downpvs - obj.thickness.median_downpvs;
            %%
            obj.thickness.epschanges = obj.thickness.bvchanges + obj.thickness.pvschanges_total;
            %%
        end

        function getdisplacement(obj)
            % GETDISPLACEMENT Calculates displacement and subtracts slow component
            obj.displacement = struct();
            % 2. Calculate and subtract slow component (using large window median filter)
            % "entire length of data" implies a very slow trend. Using 500 pts (~15-50s depending on Hz)
            obj.displacement.slow_uppvs = medfilt1(obj.idx.clean_pvsupedge_idx, 3000, 'truncate');
            obj.displacement.slow_upbv = medfilt1(obj.idx.clean_upperBVboundary, 3000, 'truncate');
            obj.displacement.slow_downbv = medfilt1(obj.idx.clean_lowerBVboundary, 3000, 'truncate');
            obj.displacement.slow_downpvs = medfilt1(obj.idx.clean_pvsdownedge_idx, 3000, 'truncate');


            % 3. Calculate changes (High-pass). Sign convention: + is OUTWARD, so
            %    the upper boundary is negated (its row index falls as it moves
            %    out) and the lower one is not.
            %    There used to be an if obj.updynamic here with the same four
            %    lines written twice, swapping only the NAMES. The negation
            %    follows the side, not the label, so the branch was pure
            %    relabelling. Removing it is what lets displacement -- the only
            %    quantity in this tree that measures motion -- be used to TEST
            %    the thickness label instead of inheriting it. See CLAUDE_LOG.md
            obj.displacement.uppvs = -(obj.idx.clean_pvsupedge_idx - obj.displacement.slow_uppvs);
            obj.displacement.downpvs = obj.idx.clean_pvsdownedge_idx - obj.displacement.slow_downpvs;
            obj.displacement.upbv = -(obj.idx.clean_upperBVboundary - obj.displacement.slow_upbv);
            obj.displacement.downbv = obj.idx.clean_lowerBVboundary - obj.displacement.slow_downbv;
        end

        function maskstack = reconstruction(obj,kymomask)
            %%
            % Thickness must match the forward crop: analyze_affine_rotate uses
            % outputHeight = floor(width), so replicate to floor(line thickness).
            crop_thickness = floor(obj.param.line_info(3,1));
            tmp.v_thr = repmat(kymomask,[1,1,crop_thickness]);
            tmp.v_thr = permute(tmp.v_thr, [3, 1, 2]); % Re-orient to [thickness, length, time]
            %%
            maskstack = analyze_affine_reverse(tmp.v_thr,obj.param.input_size,obj.param.line_info(1:2,:));
        end

        function overlay = reconstruction_overlay(obj, ch1, ch2, opts)
            %RECONSTRUCTION_OVERLAY  Burn reconstructed FWHM boundaries onto a 2-photon stack.
            %
            %   overlay = obj.RECONSTRUCTION_OVERLAY(ch1, ch2) reconstructs the BV and
            %   PVS FWHM up/down boundary lines from kymograph space back into image
            %   space (via obj.reconstruction), removes those boundary pixels from the
            %   two source channels, and returns a 4-channel hyperstack:
            %       channel 1 = ch1 with the boundary pixels blanked to 0
            %       channel 2 = ch2 with the boundary pixels blanked to 0
            %       channel 3 = reconstructed BV  boundary
            %       channel 4 = reconstructed PVS boundary
            %   The result is [ny, nx, frames, channels] — the layout io_postsavetiff wants.
            %
            %   overlay = obj.RECONSTRUCTION_OVERLAY(..., Name=Value) options:
            %       BoundaryIntensity  value written into the BV/PVS boundary channels
            %                          (default 2048; e.g. 2^16 for full 16-bit range)
            %       SavePath           folder or full '*.tif' path. Empty => do not save
            %                          (default ''). A folder gets 'fwhm_boundary.tif'.
            %       Pixel2um           pixel size in microns for the TIFF resolution tag
            %       Fps                frames per second for the TIFF finterval tag
            %
            %   Requires that fwhm()/pvsanalysis() have run so obj.mask has fields
            %   upline, downline, pvs_upline and pvs_downline. The caller sets the
            %   channel order by the order it passes ch1/ch2.
            arguments
                obj
                ch1 {mustBeNumeric}
                ch2 {mustBeNumeric}
                opts.BoundaryIntensity (1,1) double = 2048
                opts.SavePath = ''
                opts.Pixel2um (1,1) double = 1
                opts.Fps (1,1) double = 1
                opts.BlankBoundary (1,1) logical = true   % subtract boundary pixels from images for contrast
            end

            if ~isequal(size(ch1), size(ch2))
                error('line_fwhm:reconstruction_overlay:sizeMismatch', 'ch1 and ch2 must be the same size.');
            end
            reqfields = {'upline', 'downline', 'pvs_upline', 'pvs_downline'};
            missing = reqfields(~isfield(obj.mask, reqfields));
            if ~isempty(missing)
                error('line_fwhm:reconstruction_overlay:missingMask', ...
                    'obj.mask is missing field(s): %s. Run fwhm()/pvsanalysis() first.', strjoin(missing, ', '));
            end

            % 1. Reconstruct BV/PVS boundary masks (kymograph space -> image space over time).
            bv_boundary  = obj.reconstruction(obj.mask.upline     + obj.mask.downline);
            pvs_boundary = obj.reconstruction(obj.mask.pvs_upline + obj.mask.pvs_downline);

            % 2. Assemble the 4-channel stack. Optionally blank the boundary pixels
            %    out of the image channels (subtract them) to boost contrast against
            %    the coloured boundary overlay -- toggled by BlankBoundary.
            overlay = zeros([size(ch1, 1:3), 4]);
            overlay(:, :, :, 1) = ch1;
            overlay(:, :, :, 2) = ch2;

            if opts.BlankBoundary
                blank = repmat((bv_boundary + pvs_boundary) > 0, [1, 1, 1, 4]);
                overlay(blank) = 0;
            end

            overlay(:, :, :, 3) = bv_boundary  * opts.BoundaryIntensity;
            overlay(:, :, :, 4) = pvs_boundary * opts.BoundaryIntensity;

            % 3. Optionally write an ImageJ hyperstack TIFF ([ny nx frames channels]).
            savefile = char(opts.SavePath);
            if ~isempty(savefile)
                [~, ~, ext] = fileparts(savefile);
                if isempty(ext)
                    savefile = fullfile(savefile, 'fwhm_boundary.tif');
                end
                io_postsavetiff(overlay, savefile, [opts.Pixel2um, opts.Pixel2um, opts.Fps]);
            end
        end

        function save2disk(obj,name,savepath)
            line_fwhm = obj;
            save(fullfile(savepath,[name,'.mat']),'line_fwhm')
        end


    end

    methods (Static, Access = private)
        function mergedStruct = mergestruct(struct1, struct2)
            mergedStruct = struct1;
            fields2 = fieldnames(struct2);
            for i = 1:numel(fields2)
                mergedStruct.(fields2{i}) = struct2.(fields2{i});
            end
        end
    end
end

