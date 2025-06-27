classdef roi
    properties
        roimod = struct();   % struct: name -> mode ("polygon" etc.)
        vertices = struct(); % struct: name -> [x, y]
        refslice = struct(); % struct: name -> slice #
        roislice = struct(); % struct: name -> image
        stacks = struct();   % struct: name -> ROI-applied image stack
        mean = struct();
        mask = struct();
        info
        analog
    end

    methods
        function obj = roi(ref_stack, roiname, roimod, mdfExtractLoader_instance)
            arguments
                ref_stack (:,:,:) {mustBeNumeric}
                roiname (1,:) char
                roimod (1,:) char {mustBeMember(roimod, ["polygon", "rectangle"])}
                mdfExtractLoader_instance mdfExtractLoader = []
            end

            if ~isempty(mdfExtractLoader_instance)
                obj.info = mdfExtractLoader_instance.info;
                obj.analog = mdfExtractLoader_instance.analog;
            else
                obj.info = struct();
                obj.analog.data = struct();
                obj.analog.info = struct();
            end
            obj = obj.addroi(ref_stack, roiname, roimod);
        end
        
        function [obj, removed] = removeroi(obj, pattern)  
            arguments
                obj
                pattern (1,:) char
            end
        
            % Find existing ROI names                                             
            allNames = fieldnames(obj.mask);
            if isempty(allNames)
                warning('No ROIs in this object. Nothing to remove.'); %#ok<*WNTAG>
                removed = {};
                return;
            end
            % Case-insensitive substring match                                    
            idx = contains(lower(allNames), lower(pattern));
            matchNames = allNames(idx);
            % Nothing matched -- offer suggestions                                
            if isempty(matchNames)
                % Compute Levenshtein distance if available; otherwise string length diff
                try
                    dist  = cellfun(@(n) stringdistance(lower(pattern), lower(n), ...
                                                        'Method','levenshtein'), allNames);
                catch
                    dist  = abs(cellfun(@numel, allNames) - numel(pattern));  % crude fallback
                end
                [~, ord] = sort(dist);
                suggestions = strjoin(allNames(ord(1:min(3,end))), ', ');
                error('No ROI contains "%s". Did you mean: %s ?', pattern, suggestions);
            end
        
            % Remove the matched names from every ROI-indexed property            
            roiProps = {'roimod','vertices','refslice','roislice', ...
                        'stacks','mask','mean'};      % expand if you add new props
            for k = 1:numel(matchNames)
                n = matchNames{k};
                for p = 1:numel(roiProps)
                    prop = roiProps{p};
                    if isfield(obj.(prop), n)
                        obj.(prop) = rmfield(obj.(prop), n);
                    end
                end
            end
        
            removed = matchNames;  % return the list for bookkeeping, if wanted
        end


        function obj = addroi(obj, ref_stack, roiname, roimod)
                arguments
                    obj
                    ref_stack (:,:,:) {mustBeNumeric}
                    roiname (1,:) char
                    roimod (1,:) char {mustBeMember(roimod, {'rectangle','polygon','line'})}
                end
            if isfield(obj.vertices, roiname)
                fprintf('field name already exist')
                return;
            end
            vis_stack = pre_groupaverage(ref_stack,10);
            selected_roi = roiSelector(vis_stack, roimod);
            [vtx, currentsli, roimask] = selected_roi.getROI();
            obj.mask.(roiname) = roimask;
            disp(sum(roimask,'all'))
            obj.roimod.(roiname) = roimod;
            obj.vertices.(roiname) = vtx;
            obj.roislice.(roiname) = vis_stack(:,:,currentsli);
            obj.refslice.(roiname) = round(currentsli*10);
        end

        function obj = modifyroi(obj, ref_stack, roiname)
            vis_stack = pre_groupaverage(ref_stack,10);
            selected_roi = roiSelector(vis_stack, obj.roimod.(roiname),obj.vertices.(roiname));
            [vtx, currentsli,roimask] = selected_roi.getROI();
            disp(sum(obj.mask.(roiname),'all'))
            obj.mask.(roiname) = roimask;
            disp(sum(roimask,'all'))
            obj.vertices.(roiname) = vtx;
            obj.roislice.(roiname) = vis_stack(:,:,currentsli);
            obj.refslice.(roiname) = round(currentsli*10);
        end

        function roistack = addstack(obj, stack, roiname)
            roistack = roi_applyvertices(stack, obj.vertices.(roiname));
        end

        function overlay = showroi(obj, roiname)
            overlay = roi_summaryimg(obj.roislice.(roiname), obj.vertices.(roiname));
        end

        function state = showstack(obj, stackfieldname, roiname)
            state = util_checkstack(obj.stacks.(roiname).(stackfieldname));
        end

        function obj = copyroi(obj, original_name, new_name)
            % Duplicate all fields for a given ROI name
            if ~isfield(obj.vertices, original_name)
                error("ROI '%s' does not exist.", original_name);
            end
            % Copy each ROI-related field
            obj.roimod.(new_name)   = obj.roimod.(original_name);
            obj.vertices.(new_name) = obj.vertices.(original_name);
            obj.refslice.(new_name) = obj.refslice.(original_name);
            obj.roislice.(new_name) = obj.roislice.(original_name);
            obj.mask.(new_name) = obj.mask.(original_name);
        end
      
    end % end of public method  



properties (Constant)
            clist = struct(...
            'white', [1,1,1],...
            'black', [0,0,0],...
            'cyan', [0,1,1],...
            'magenta', [1,0,1],...
            'yellow', [1,1,0],...
            'red', [1,0,0],...
            'green', [0,1,0],...
            'blue', [0,0,1],...
            'orange', [1,0.5,0],...
            'lightgreen', [0.4,0.8,0.4],...
            'darkgreen', [0,0.5,0]);
    end
end
