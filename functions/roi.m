classdef roi
    %ROI Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mod
        vertices
        refslice
        roislice 
        stacks = struct();
        data = struct();
    end
    
    methods
        function obj = roi(ref_stack,mod)
            obj.mod = mod;
            vis_stack = pre_groupaverage(ref_stack,10);
            [obj.vertices,currentsli] = roi_rectangle_polygon(vis_stack,mod);
            obj.roislice = vis_stack(:,:,currentsli);
            obj.refslice = round(currentsli*10);
        end
% modifyroi function is intended to fine adjustment using multiple color
% channel as its hard to 
        function obj = modifyroi(obj,ref_stack) 
            vis_stack = pre_groupaverage(ref_stack,10);
            [obj.vertices,currentsli] = roi_rectangle_polygon(vis_stack,obj.mod,obj.vertices);
            obj.roislice = vis_stack(:,:,currentsli);
            obj.refslice = round(currentsli*10);
            obj.stacks = struct();
            obj.data = struct();
        end
        % applyvertices simply
        function roistack = addstack(obj,stack)
            roistack = roi_applyvertices(stack,obj.vertices);
        end
        %
        function overlay = showroi(obj)
            overlay = roi_summaryimg(obj.roislice,obj.vertices);
        end
        
        function state = showstack(obj,stackfieldname)
            state = util_checkstack(getfield(obj.stacks,stackfieldname));
        end

        function obj = radonthresholding(obj,stackfieldname)
            stack = getfield(obj.stacks,stackfieldname);
            stack(isnan(stack))=0;
            filteredStack = gpuArray(stack);
            for i = 1:size(filteredStack, 3)
                sli = medfilt2(filteredStack(:, :, i), [3 3]);  % Apply 3x3 median filter (adjust size as needed)
                filteredStack(:, :, i) = imgaussfilt(sli,1);
            end
            filteredStack = gather(filteredStack);
            % filteredStack = pre_thresholding(filteredStack);
            irtd_fieldname = ['irtd_',stackfieldname];
            [obj.stacks.(irtd_fieldname),obj.stacks.(tirs_fieldname)] = analyze_radon(filteredStack);
            end
        end

end

