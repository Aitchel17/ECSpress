classdef roi
    % roi class can constructed with reference_stack
    
    properties
        roimod
               
        vertices % vertices for this roi
        refslice % reference slice number when set a roi
        roislice % slice used for roi creation
        stacks = struct(); % stacks of roi
        data = struct(); % data extracted from roi stack

        info % roi's metadata
        analog %
    end
    
    methods
        function obj = roi(reference_stack,roimod,mdfExtractLoader_instance)
            arguments
                reference_stack (:,:,:) {mustBeNumeric}
                roimod (1,:) char {mustBeMember(roimod, ["polygon", "rectangle"])}
                mdfExtractLoader_instance mdfExtractLoader = []
            end
            % if the user is using mdfExtractLoader,
            %  use predifinedmetadata structure
            if ~isempty(mdfExtractLoader_instance)
                obj.info = mdfExtractLoader_instance.info;
                obj.analog = mdfExtractLoader_instance.analog;
            else
                obj.info = struct();
                obj.analog.data = struct();
                obj.analog.info = struct();
            end
            obj.roimod = roimod;
            vis_stack = pre_groupaverage(reference_stack,10);
            [obj.vertices,currentsli] = roi_rectangle_polygon(vis_stack,roimod);
            obj.roislice = vis_stack(:,:,currentsli);
            obj.refslice = round(currentsli*10);
        end


% modifyroi function is intended to fine adjustment using multiple color
% channel as its hard to 
        function obj = modifyroi(obj,ref_stack) 
            vis_stack = pre_groupaverage(ref_stack,10);
            [obj.vertices,currentsli] = roi_rectangle_polygon(vis_stack,obj.roimod,obj.vertices);
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


    end
end

