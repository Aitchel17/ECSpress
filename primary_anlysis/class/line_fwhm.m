classdef line_fwhm < handle
    %LINE_FWHM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        rotatecrop = struct()
        kymograph = struct()
        mask
        idx
        t_axis
        param = struct('bv_thr',-1,'bv_offset',-1,'csf_thr',-1,'csf_offset',-1,...
            'input_size',[-1,-1,-1], 'line_info',[-1,-1;-1,-1;-1,-999]); % 3x2 double, {x1,y1; x2,y2; linewidth,-999}
    end
    
    
    methods
        function obj = line_fwhm(line_info)
            % Constructor function get roi info, ex: (x1,y1;x2,y2;,thickness,-999)
            arguments
                line_info (3,2)
            end
            obj.param.line_info = line_info;
        end

        function addkymograph(obj,stack_name,stack)
            arguments
                obj
                stack_name (1,1) string {mustBeMember(stack_name,{'lumen','wall', 'pvs', 'outside'})}
                stack
            end
            obj.param.input_size = size(stack);
            name_rotatecrop = strcat("rc_",stack_name);
            name_kymograph = strcat("kgph_",stack_name);
            name_normkymograph = strcat("normkgph_",stack_name);
            rotatecroped_stack =analyze_affine_rotate(stack,obj.param.line_info(1:2,:), obj.param.line_info(3,1));
            rawkymograph = squeeze(sum(rotatecroped_stack, 1));
            obj.rotatecrop.(name_rotatecrop) = rotatecroped_stack;
            obj.kymograph.(name_kymograph) = rawkymograph;
            obj.kymograph.(name_normkymograph) = (rawkymograph-min(rawkymograph,[],1))./max(rawkymograph,[],1); % just for reference to see intermediate step, no offset
        end
        
        function kymograph_afterprocess(obj,stack_name,window_size)
         arguments
                obj
                stack_name (1,1) string {mustBeMember(stack_name,{'lumen','wall', 'pvs', 'outside'})}
                window_size   (1,2) double = [3 2] % default halfmax
         end
            name_kymograph = strcat("kgph_",stack_name);
            name_processed_kymograph = strcat(name_kymograph,'_processed');
            processed_kymograph = medfilt2(obj.kymograph.(name_kymograph),window_size);
            obj.kymograph.(name_processed_kymograph) = processed_kymograph;

        end


        function fwhm(obj,stack_name,threshold,offset)
             arguments
                obj
                stack_name (1,1) string {mustBeMember(stack_name,{'lumen','wall', 'pvs', 'outside'})}
                threshold   (1,1) double = 0.5 % default halfmax
                offset      (1,1) double = 10 % default 10% offset
             end
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here     
            name_kymograph = strcat("kgph_",stack_name,"_processed");
            [tmp.idx, tmp.kgph_mask] = analyze_fwhm(obj.kymograph.(name_kymograph),threshold,offset);
            obj.param.bv_thr = threshold;
            obj.param.bv_offset = offset;
            obj.idx = obj.mergestruct(obj.idx, tmp.idx);
            obj.mask = obj.mergestruct(obj.mask, tmp.kgph_mask);
        end

        function csfanalysis(obj,csf_stack,threshold,offset)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            obj.rotatecrop.rc_csf =analyze_affine_rotate(csf_stack,obj.param.line_info(1:2,:), obj.param.line_info(3,1));
            obj.kymograph.kgph_csf = squeeze(sum(obj.rotatecrop.rc_csf,1));
            obj.kymograph.kgph_normcsf = (obj.kymograph.kgph_csf-min(obj.kymograph.kgph_csf,[],1))./max(obj.kymograph.kgph_csf,[],1);

             [tmp.idx, tmp.kgph_mask] = analyze_csfoutter(obj.kymograph.kgph_csf, obj.idx.upperboundary, obj.idx.lowerboundary, threshold,offset);

             obj.param.csf_thr = threshold;
             obj.param.csf_offset = offset;
             obj.idx = obj.mergestruct(obj.idx, tmp.idx);
             obj.mask = obj.mergestruct(obj.mask, tmp.kgph_mask);
        end

        function maskstack = reconstruction(obj,kymomask)
            tmp.v_thr = repmat(kymomask,[1,1,size(obj.rotatecrop.rc_bv,1)]);
            tmp.v_thr = permute(tmp.v_thr,[3,1,2]);
            maskstack = analyze_affine_reverse(tmp.v_thr,obj.param.input_size,obj.param.line_info(1:2,:));
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

