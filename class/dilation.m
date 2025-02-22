classdef dilation < roi
    %EPS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = dilation(stackname,roitype)
            obj@roi(stackname,roitype)
        end
        

        function obj = radonthresholding(obj,stackfieldname)
            disp([stackfieldname ' radonthresholding started'])
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
            tirs_fieldname = ['tirs_',stackfieldname];
            [obj.stacks.(irtd_fieldname),obj.stacks.(tirs_fieldname)] = analyze_radon(filteredStack); % do tirs
        end

    end
end

