classdef roi_handle < handle
    properties
        ROIs = struct('Label',{}, 'Mode',{}, 'Vertices',{}, ...
            'Mask',{}, 'ImageSize',{}, ...
            'RefSlice',{}, 'ROISlice',{}, ...
            'Created',{}, 'Modified',{});
        savepath (1,:) char
        isloaded logical
    end

    methods
        function obj = roi_handle(savepath)
            if  isfile(savepath)
                disp('previous roilist detected loading the roilist')
                loadstruct = load(savepath);
                obj = loadstruct.roilist;
                obj.isloaded = true;
            else
                obj.savepath = savepath;
                obj.isloaded = false;
            end
        end

        function addormodifyroi(obj, ref_stack, label, roimode)
            arguments
                obj
                ref_stack (:,:,:) {mustBeNumeric}
                label (1,:) char
                roimode (1,:) char {mustBeMember(roimode, {'rectangle','polygon','line'})}
            end
            try
                obj.addroi(ref_stack, label, roimode);
            catch
                disp(['ROI "' label '" already exists. Modifying...']);
                obj.modifyroi(ref_stack, label);
            end
        end

        function addroi(obj, ref_stack, label, roimode)
            arguments
                obj
                ref_stack (:,:,:) {mustBeNumeric}
                label   (1,:) char
                roimode (1,:) char {mustBeMember(roimode, {'rectangle','polygon','line'})}
            end
            % 라벨 중복 금지
            if any(strcmp({obj.ROIs.Label}, label))
                error('ROI label "%s" already exists.', label);
            end

            vis_stack = pre_groupaverage(ref_stack, 10);
            selector  = roiSelector(vis_stack, roimode);
            [vtx, sli, roimask] = selector.getROI();

            nowt = datetime('now');
            newroi.Label     = label;
            newroi.Mode      = roimode;
            newroi.Vertices  = vtx;
            newroi.Mask      = roimask;
            newroi.ImageSize = size(roimask);
            newroi.RefSlice  = floor(sli*10);
            newroi.ROISlice = [];
            newroi.Created   = nowt;
            newroi.Modified  = nowt;
            obj.ROIs(end+1) = newroi;
        end
        function addrgb(obj, rgb_img, label, roimode)
            arguments
                obj
                rgb_img (:,:,:) {mustBeNumeric}
                label   (1,:) char
                roimode (1,:) char {mustBeMember(roimode, {'rectangle','polygon','line'})}
            end
            if any(strcmp({obj.ROIs.Label}, label))
                error('ROI label "%s" already exists.', label);
            end

            selector  = roiSelector(rgb_img, roimode);
            [vtx, sli, roimask] = selector.getROI();

            nowt = datetime('now');
            newroi.Label     = label;
            newroi.Mode      = roimode;
            newroi.Vertices  = vtx;
            newroi.Mask      = roimask;
            newroi.ImageSize = size(roimask);
            newroi.RefSlice  = floor(sli*10);
            newroi.ROISlice = [];
            newroi.Created   = nowt;
            newroi.Modified  = nowt;
            obj.ROIs(end+1) = newroi;
        end

        function addimgchannel(obj, imgstack, labels)
            idx_list = [];
            for i = 1:length(labels)
                i = obj.findLabel(labels(i));
                idx_list = [idx_list, i];
            end

            if ndims(imgstack) == 4
                roistruct = obj.ROIs;
                for i = 1:length(idx_list)
                    sli_start = obj.ROIs(idx_list(i)).RefSlice -5;
                    sli_end = obj.ROIs(idx_list(i)).RefSlice + 5;
                    if sli_start < 1
                        sli_end = sli_end - sli_start+1;
                        sli_start = 1;
                    end
                    if sli_end > size(imgstack,3)
                        sli_start = sli_start - (size(imgstack,3)-sli_end);
                        sli_end = size(imgstack,3);
                    end
                    averaged_stack = mean(imgstack(:,:,sli_start:sli_end,:),3);
                    averaged_stack = squeeze(averaged_stack);
                    obj.ROIs(idx_list(i)).ROISlice = averaged_stack;
                    %% debug
                    % rgb_image = zeros([size(averaged_stack,1,2),3]);
                    % rgb_image(:,:,1:2) = averaged_stack;
                    % rgb_image = mat2gray(rgb_image);
                    % figure()
                    % imshow(rgb_image)
                    % axis image
                    % colormap('gray')
                    %%
                end
            end
        end

        function modifyroi(obj, ref_stack, label)
            arguments
                obj
                ref_stack (:,:,:) {mustBeNumeric}
                label   (1,:) char
            end
            i = obj.findLabel(label);
            vis_stack = pre_groupaverage(ref_stack, 10);
            selector  = roiSelector(vis_stack, obj.ROIs(i).Mode, obj.ROIs(i).Vertices);
            [vtx, sli, roimask] = selector.getROI();
            obj.ROIs(i).Vertices  = vtx;
            obj.ROIs(i).Mask      = roimask;
            obj.ROIs(i).ImageSize = size(roimask);
            obj.ROIs(i).RefSlice  = round(sli*10);
            obj.ROIs(i).ROISlice  = vis_stack(:,:,sli);
            obj.ROIs(i).Modified  = datetime('now');
        end

        function modifyrgb(obj, ref_stack, label)
            arguments
                obj
                ref_stack (:,:,:) {mustBeNumeric}
                label   (1,:) char
            end
            i = obj.findLabel(label);
            selector  = roiSelector(ref_stack, obj.ROIs(i).Mode, obj.ROIs(i).Vertices);
            [vtx, sli, roimask] = selector.getROI();
            obj.ROIs(i).Vertices  = vtx;
            obj.ROIs(i).Mask      = roimask;
            obj.ROIs(i).ImageSize = size(roimask);
            obj.ROIs(i).RefSlice  = round(sli*10);
            obj.ROIs(i).ROISlice  = ref_stack(:,:,sli);
            obj.ROIs(i).Modified  = datetime('now');
        end

        function copyroi(obj, original_label, new_label)
            arguments
                obj
                original_label (1,:) char
                new_label      (1,:) char
            end
            % 원본 존재/신규 라벨 중복 체크
            i = obj.findLabel(original_label);
            if any(strcmp({obj.ROIs.Label}, new_label))
                error('ROI label "%s" already exists.', new_label);
            end

            nowt = datetime('now');
            E = obj.ROIs(i);
            E.Label    = new_label;
            E.Created  = nowt;
            E.Modified = nowt;
            % (ROISlice/Mask/Vertices 등은 E에 이미 복사돼 있음)
            obj.ROIs(end+1) = E;
        end

        function removedLabel = removeroi(obj, label)
            arguments
                obj
                label (1,:) char
            end
            i = obj.findLabel(label);
            removedLabel = obj.ROIs(i).Label;
            obj.ROIs(i) = [];
        end



        function roistack = applyvertices(obj, stack, label)
            arguments
                obj
                stack (:,:,:) {mustBeNumeric}
                label (1,:) char
            end
            i = obj.findLabel(label);
            roistack = roi_applyvertices(stack, obj.ROIs(i).Vertices);
        end

        function vertices = getvertices(obj, label)
            i = obj.findLabel(label);
            vertices = obj.ROIs(i).Vertices;
        end

        function vertices = getmask(obj, label)
            i = obj.findLabel(label);
            vertices = obj.ROIs(i).Mask;
        end

        function labels = list(obj)
            labels = string({obj.ROIs.Label});
        end

        function save2disk(obj)
            roilist = obj;
            save(obj.savepath,'roilist')
        end

    end

    methods (Access=?make_fig)
        function i = findLabel(obj, label)
            if isempty(obj.ROIs)
                error('No ROI stored.');
            end
            i = find(strcmp({obj.ROIs.Label}, label), 1);
            if isempty(i)
                error('ROI "%s" not found.', label);
            end
        end
    end
end
