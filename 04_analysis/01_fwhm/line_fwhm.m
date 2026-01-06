classdef line_fwhm < handle
    %LINE_FWHM Summary of this class goes here
    %   Detailed explanation goes here

    properties
        rotatecrop = struct()
        kymograph = struct()
        thickness = struct()
        displacement = struct()
        updynamic
        mask
        idx
        t_axis
        param = struct('input_size',[-1,-1,-1], 'line_info',[-1,-1;-1,-1;-1,-999]); % 3x2 double, {x1,y1; x2,y2; linewidth,-999}
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
            name_kymograph = strcat("kgph_",stack_name);
            name_processed_kymograph = strcat(name_kymograph,'_processed');
            processed_kymograph = medfilt2(obj.kymograph.(name_kymograph),window_size);
            prctile10 = prctile(processed_kymograph,10,1);
            prctile90 = prctile(processed_kymograph,90,1);
            processed_kymograph = min(processed_kymograph, prctile90);
            processed_kymograph = max(processed_kymograph, prctile10);
            processed_kymograph = processed_kymograph - min(processed_kymograph,[],1);
            processed_kymograph = processed_kymograph./max(processed_kymograph,[],1);
            obj.kymograph.(name_processed_kymograph) = processed_kymograph;
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

        function clean_outlier(obj,overwrite)
            arguments
                obj
                overwrite = false
            end
            obj.idx = clean_outlier(obj.idx, overwrite);
        end

        function getdiameter(obj)
            obj.thickness = struct();
            obj.thickness.bv = obj.idx.clean_lowerBVboundary - obj.idx.clean_upperBVboundary;
            uppvs_thickness = obj.idx.clean_upperBVboundary - obj.idx.clean_pvsupedge_idx;
            downpvs_thickness = obj.idx.clean_pvsdownedge_idx - obj.idx.clean_lowerBVboundary;
            obj.thickness.totalpvs = uppvs_thickness + downpvs_thickness;
            difference_pvs = uppvs_thickness - downpvs_thickness;
            difference_pvs = medfilt1(difference_pvs,11); % smoothing the pvs thickness difference
            % caclulate area under the curve of difference_pvs
            difference_pvs = trapz(difference_pvs);
            if difference_pvs > 0
                disp('up pvs is dynamic')
                obj.updynamic = true;
                obj.thickness.dynamic_pvs = uppvs_thickness;
                obj.thickness.static_pvs = downpvs_thickness;
            else
                disp('up pvs is static')
                obj.updynamic = false;
                obj.thickness.dynamic_pvs = downpvs_thickness;
                obj.thickness.static_pvs = uppvs_thickness;
            end

            %%
            obj.thickness.median_totalpvs = median(obj.thickness.totalpvs);
            obj.thickness.median_staticpvs = median(obj.thickness.static_pvs);
            obj.thickness.median_dynamicpvs = median(obj.thickness.dynamic_pvs);
            %%
            obj.thickness.median_bv = median(obj.thickness.bv);
            obj.thickness.bvchanges = obj.thickness.bv - obj.thickness.median_bv;
            obj.thickness.pvschanges_total = obj.thickness.totalpvs - obj.thickness.median_totalpvs;
            obj.thickness.pvschanges_dynamic = obj.thickness.dynamic_pvs - obj.thickness.median_dynamicpvs;
            obj.thickness.pvschanges_static = obj.thickness.static_pvs - obj.thickness.median_staticpvs;
            %%
            obj.thickness.ecschanges_residual = obj.thickness.bvchanges + obj.thickness.pvschanges_total;
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


            % 3. Calculate changes (High-pass)
            if obj.updynamic
                obj.displacement.dynamicpvs = -(obj.idx.clean_pvsupedge_idx - obj.displacement.slow_uppvs);
                obj.displacement.staticpvs = obj.idx.clean_pvsdownedge_idx - obj.displacement.slow_downpvs;
                obj.displacement.dynamicbv = -(obj.idx.upperBVboundary - obj.displacement.slow_upbv);
                obj.displacement.staticbv = obj.idx.lowerBVboundary - obj.displacement.slow_downbv;
            else
                obj.displacement.staticpvs = -(obj.idx.clean_pvsupedge_idx - obj.displacement.slow_uppvs);
                obj.displacement.dynamicpvs = obj.idx.clean_pvsdownedge_idx - obj.displacement.slow_downpvs;
                obj.displacement.staticbv = -(obj.idx.upperBVboundary - obj.displacement.slow_upbv);
                obj.displacement.dynamicbv = obj.idx.lowerBVboundary - obj.displacement.slow_downbv;
            end
        end

        function maskstack = reconstruction(obj,kymomask)
            tmp.v_thr = repmat(kymomask,[1,1,size(obj.rotatecrop.rc_bv,1)]);
            tmp.v_thr = permute(tmp.v_thr,[3,1,2]);
            maskstack = analyze_affine_reverse(tmp.v_thr,obj.param.input_size,obj.param.line_info(1:2,:));
        end

        function save2disk(obj,savepath)
            line_fwhm = obj;
            save(fullfile(savepath,'line_fwhm.mat'),'line_fwhm')
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

