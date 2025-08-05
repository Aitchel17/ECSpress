classdef line_fwhm
    %LINE_FWHM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        rotatecrop = struct()
        kymograph = struct('kgph_bv',[],'kgph_csf',[],'kgph_normcsf',[],'kgph_normbv',[])
        mask
        idx
        t_axis
        param = struct('bv_thr',-1,'bv_offset',-1,'csf_thr',-1,'csf_offset',-1,...
            'input_size',[-1,-1,-1], 'line_info',[-1,-1;-1,-1;-1,-999]); % 3x2 double, {x1,y1; x2,y2; linewidth,-999}
    end
    
    methods
        function obj = line_fwhm(vessel_stack,csf_stack,line_info)
            obj.param.line_info = line_info;
            obj.param.input_size = size(vessel_stack);
            obj.rotatecrop.rc_bv =analyze_affine_rotate(vessel_stack,line_info(1:2,:), line_info(3,1));
            obj.rotatecrop.rc_csf =analyze_affine_rotate(csf_stack,line_info(1:2,:), line_info(3,1));
            obj.kymograph.kgph_bv = squeeze(sum(obj.rotatecrop.rc_bv,1));
            obj.kymograph.kgph_csf = squeeze(sum(obj.rotatecrop.rc_csf,1));
            obj.kymograph.kgph_normbv = (obj.kymograph.kgph_bv-min(obj.kymograph.kgph_bv,[],1))./max(obj.kymograph.kgph_bv,[],1);
            obj.kymograph.kgph_normcsf = (obj.kymograph.kgph_csf-min(obj.kymograph.kgph_csf,[],1))./max(obj.kymograph.kgph_csf,[],1);
        end
        
        function obj = bvanalysis(obj,threshold,offset)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
             [tmp.idx, tmp.kgph_mask] = analyze_fwhm(obj.kymograph.kgph_bv,threshold,offset);
             obj.param.bv_thr = threshold;
             obj.param.bv_offset = offset;

             obj.idx = obj.mergestruct(obj.idx, tmp.idx);
             obj.mask = obj.mergestruct(obj.mask, tmp.kgph_mask);
        end

        function obj = csfanalysis(obj,threshold,offset)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
             [tmp.idx, tmp.kgph_mask] = analyze_csfoutter(obj.kymograph.kgph_csf,obj.idx.bv_upperboundary,obj.idx.bv_lowerboundary,threshold,offset);
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

