classdef line_fwhm < handle
    %LINE_FWHM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties 
        rotatecrop = struct()
        kymograph = struct()
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
            [tmp.idx, tmp.kgph_mask,tmp.param] = analyze_fwhm(obj.kymograph.(name_kymograph));
            % obj.param.bv_thr = threshold;
            % obj.param.bv_offset = offset;
            obj.idx = obj.mergestruct(obj.idx, tmp.idx);
            obj.mask = obj.mergestruct(obj.mask, tmp.kgph_mask);
        end

        function pvsanalysis(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
             [tmp.idx, tmp.kgph_mask] = analyze_csfoutter(obj.kymograph.kgph_pvs_processed, obj.idx.upperboundary, obj.idx.lowerboundary);
             obj.idx = obj.mergestruct(obj.idx, tmp.idx);
             obj.mask = obj.mergestruct(obj.mask, tmp.kgph_mask);
        end

        function clean_outlier(obj,overwrite)
            idx_dataset = obj.idx;
            idx_names = fieldnames(idx_dataset);
            idx_names = string(idx_names');
            processed_nameidx = startsWith(idx_names,'clean_');
            processed_names = idx_names(processed_nameidx);
            processed_names = erase(processed_names,'clean_');
            original_names = idx_names(~processed_nameidx);
            if overwrite
                unprocessed_names = original_names;
            else
                match_idx = ismember(processed_names,original_names);
                unprocess_logic = true(size(original_names));
                unprocess_logic(match_idx) = false;
                unprocessed_names = original_names(unprocess_logic);
            end

            %%
            for idx_name = unprocessed_names
                if isempty(idx_name)
                    disp('all processed')
                else
                    disp(idx_name{:})
                    idx_data = idx_dataset.(idx_name{:});
                    ref_median = medfilt1(idx_data,31);
                    diff_data = abs(idx_data - ref_median);
                    above_noise = diff_data > 4;
                    % above data already contains the positional information
                    % make edge case at the starting point and the end point
                    % Can not interpolate at the start and end
                    % if it differ, just put median value
                    % idx matching is tricky.. the interpolation 
                    edge_thresholded = above_noise;
                    edge_thresholded(1) = 0; % ensure that transit start from rising edge, if start is 1 falling edge is longer
                    edge_thresholded(end) = 0; % ensure that transit end with falling edge, if end is 1 rising edge would longer
                    transit_thr = diff(edge_thresholded);
                    %
                    if sum(transit_thr) == 0
                        rise_logic = transit_thr == 1; 
                        rise_logic = [rise_logic,false];
                        fall_logic = transit_thr == -1;
                        fall_logic = [false,fall_logic];
                    else
                        disp('number of rising and falling is different')
                    end
                    interp_data = (idx_data(rise_logic)+idx_data(fall_logic))/2;
                    %
                    run = [1, abs(transit_thr)];
                    run = cumsum(run);  % make transition point to have unique value
                    unique_run = run(edge_thresholded==1); % extract values at the threshold
                    unique_run = unique(unique_run);   % get unique values
                    [~,runidx] = ismember(run,unique_run); % find position of unique run if run matches unique_run
                    
                    idxdata_interp_idx = runidx>0;
                    interp_idx = runidx(idxdata_interp_idx);
                    
                    corrected_data = idx_data;
                    %
                    corrected_data(idxdata_interp_idx) = interp_data(interp_idx);  
                    idx_dataset.(['clean_' idx_name{:}]) = corrected_data;
                    obj.idx = idx_dataset;
                end
                
            end
            %%
          
            
            %% debuging code
            % tfig = figure();
            % tax = axes();
            % %%
            % figure(tfig)
            % cla(tax)
            % hold(tax,'on')
            % plot(rise_logic+mean(ref_median),'r*')
            % plot(fall_logic+mean(ref_median),'go')
            % plot(edge_thresholded+mean(ref_median),'ms')
            % %%
            % plot(line_data,'g')
            % plot(ref_median,'r')
            % plot(corrected_data,'k')
            %%

            %%

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

