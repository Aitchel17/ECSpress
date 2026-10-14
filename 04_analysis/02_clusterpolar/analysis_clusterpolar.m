classdef analysis_clusterpolar < handle
%ANALYSIS_CLUSTERPOLAR  Kymograph clusters, three diameter states, their contours and polar profiles.
%   Caller: analysis_main, analysis_main_grabne, analysis_phantom, fig_review_session,
%           polarcluster_integration, centralize_strip.centralize_strippolar, ECSSession
%
%   pc = analysis_clusterpolar();
%   pc.cluster(kgph_lumen, stack_bv, stack_pvs, upper_bv, lower_bv);   % seconds to minutes
%   pc.pickstate();                              % two picks by hand, intermediate by rule
%   pc.drawcontour(roilist);                     % six selection polygons by hand
%   pc.transform(roilist);                       % contours -> polar profiles
%   pc.save2disk(primary_dir);                   % polarcluster.mat, the object itself
%   pc = analysis_clusterpolar.from_struct(s);   % a polarcluster.mat written before the class

    properties
        param                       % struct        cluster_num (count)
        filt                        % struct        min_duration (frames)
        kgph_lumen_columnidx        % 1 x T int     frame order that makes each cluster contiguous
        clusterboundary             % K x 2 int     first and last position of each cluster in that order
        clust_med_bv                % H x W x K     per-cluster median image, lumen channel
        clust_med_csf               % H x W x K     per-cluster median image, PVS channel
        cluster_summary             % 1 x K struct  id, mean_diameter, std_diameter (px), duration (frames), mean_lowchange, mean_upchange
        filtered_clusteridx         % 1 x k int     clusters longer than filt.min_duration
        constricted_cluster_idx     % 1 x 1 int     picked by hand
        dilated_cluster_idx         % 1 x 1 int     picked by hand
        intermediate_cluster_idx    % 1 x 1 int     chosen by chooseintermediate
        polar_theta                 % 1 x 360 rad   bin centres, shared by every profile
        polar_profiles              % 1 x 6 struct  roi_label, binned_r (1 x 360 px, NaN = empty bin); c/d/i x bv/pvs
        polar_binstart_deg          % 1 x 1 deg     added to polar_theta to get the image angle
        pax_angle                   % 1 x 1 deg     PAX line, display only
    end

    properties (Dependent)
        constricted_medianimg       % H x W         PVS-channel median of the constricted cluster
        dilated_medianimg           % H x W         PVS-channel median of the dilated cluster
    end

    methods
        function obj = analysis_clusterpolar()
            obj.param.cluster_num = 20;
            obj.filt.min_duration = 10;
        end

        function cluster(obj, kgph_lumen, stack_bv, stack_pvs, upper_bv, lower_bv)
        % Cluster the lumen kymograph, median each cluster's frames in both stacks, summarise and filter
        %
        % IN   kgph_lumen  P x T float   processed lumen kymograph, one column per frame
        %      stack_bv    H x W x T     lumen-channel stack, the caller resolves which channel
        %      stack_pvs   H x W x T     PVS-channel stack
        %      upper_bv    1 x T px      upper lumen wall per frame, for the diameter summary
        %      lower_bv    1 x T px      lower lumen wall per frame
            [obj.kgph_lumen_columnidx, obj.clusterboundary] = analysis_cluster_kymograph( ...
                kgph_lumen, obj.param.cluster_num);
            sorted_bv = stack_bv(:, :, obj.kgph_lumen_columnidx);
            sorted_csf = stack_pvs(:, :, obj.kgph_lumen_columnidx);

            n_cluster = obj.param.cluster_num;
            obj.clust_med_bv = zeros([size(sorted_bv, 1, 2), n_cluster]);
            obj.clust_med_csf = zeros([size(sorted_csf, 1, 2), n_cluster]);
            for cluster_idx = 1:n_cluster
                span = obj.clusterboundary(cluster_idx, 1):obj.clusterboundary(cluster_idx, 2);
                obj.clust_med_bv(:, :, cluster_idx) = median(sorted_bv(:, :, span), 3);
                obj.clust_med_csf(:, :, cluster_idx) = median(sorted_csf(:, :, span), 3);
            end

            obj.summarize(upper_bv, lower_bv);
            is_long = [obj.cluster_summary.duration] > obj.filt.min_duration;
            obj.filtered_clusteridx = [obj.cluster_summary(is_long).id];
        end

        function pickstate(obj)
        % Constricted and dilated by hand on the filtered PVS medians; intermediate by rule
            filtered_csf = obj.clust_med_csf(:, :, obj.filtered_clusteridx);
            constricted_pos = util_checkstack(filtered_csf, 'Select the constricted image.');
            dilated_pos = util_checkstack(filtered_csf, 'Select the dilated image.');
            obj.constricted_cluster_idx = obj.filtered_clusteridx(constricted_pos);
            obj.dilated_cluster_idx = obj.filtered_clusteridx(dilated_pos);
            obj.chooseintermediate();
        end

        function img = stateimage(obj, channel, state)
        % One state's cluster median, straight out of clust_med_bv / clust_med_csf
            arguments
                obj
                channel (1,1) string {mustBeMember(channel, ["bv", "pvs"])}
                state   (1,1) string {mustBeMember(state, ["constricted", "intermediate", "dilated"])}
            end
            cluster_idx = obj.(char(state + "_cluster_idx"));
            if channel == "bv"
                img = obj.clust_med_bv(:, :, cluster_idx);
            else
                img = obj.clust_med_csf(:, :, cluster_idx);
            end
        end

        function img = get.constricted_medianimg(obj)
            img = obj.stateimage("pvs", "constricted");
        end

        function img = get.dilated_medianimg(obj)
            img = obj.stateimage("pvs", "dilated");
        end

        function drawcontour(obj, roilist)
        % One selection polygon per state and channel, kept in roilist; dilated opens before
        % intermediate so an intermediate with no polygon of its own starts from it
            state_names = {'constricted', 'intermediate', 'dilated'};
            draw_order  = [1 3 2];   % 1 x 3 int   state indices in the order the windows open
            ch_names    = {'bv', 'pvs'};

            for state_idx = draw_order
                for ch_idx = 1:2
                    contour_name = sprintf('%s%s_contour', state_names{state_idx}, upper(ch_names{ch_idx}));
                    disp_msg     = sprintf('%s %s', state_names{state_idx}, ch_names{ch_idx});

                    img = obj.stateimage(ch_names{ch_idx}, state_names{state_idx});
                    % Window on the tissue, not on the zeros; see CLAUDE_LOG.md "3% of the greys"
                    disp_lim = prctile(img(img > 0), [1 99]);

                    % 1. Auto-detect edges
                    disp(['Auto-detecting edges for: ' disp_msg '...']);
                    try
                        img_smooth = imgaussfilt(double(img), 3);
                        edges = edge(img_smooth, 'log', 0.005);
                    catch ME
                        disp(['Error in auto-detection: ' ME.message]);
                        edges = false(size(img));
                    end

                    % 2. The polygon to start from: this label's own, else the dilated one
                    existing_pos = [];
                    roi_idx = [];
                    if ~isempty(roilist) && ~isempty(roilist.ROIs)
                        roi_idx = find(strcmp({roilist.ROIs.Label}, contour_name), 1);
                        if ~isempty(roi_idx)
                            existing_pos = roilist.ROIs(roi_idx).Vertices;
                            disp(['Found existing ROI in roilist: ' contour_name]);
                        end
                    end
                    is_unseeded_intermediate = isempty(roi_idx) && strcmp(state_names{state_idx}, 'intermediate');
                    if is_unseeded_intermediate && ~isempty(roilist) && ~isempty(roilist.ROIs)
                        seed_name = sprintf('dilated%s_contour', upper(ch_names{ch_idx}));
                        seed_idx = find(strcmp({roilist.ROIs.Label}, seed_name), 1);
                        if ~isempty(seed_idx)
                            existing_pos = roilist.ROIs(seed_idx).Vertices;
                            disp(['Starting from ' seed_name '; shrink it onto the intermediate wall']);
                        end
                    end

                    % 3. Figure: image, detected edges in red
                    hF = figure('Name', ['Select Edges: ' disp_msg], 'NumberTitle', 'off');
                    hIm = imshow(img, disp_lim);
                    uicontrol(hF, 'Style', 'pushbutton', 'String', 'Contrast', ...
                        'Units', 'normalized', 'Position', [0.005 0.005 0.09 0.05], ...
                        'Callback', @(~,~) imcontrast(hIm));   % its own window, so drawpolygon keeps the focus
                    hold on;
                    [Ay, Ax] = find(edges);
                    if ~isempty(Ax)
                        plot(Ax, Ay, 'r.', 'MarkerSize', 2);
                    else
                        text(10, 10, 'No edges detected', 'Color', 'y');
                    end
                    title({['Step 1: Detected Edges for ' disp_msg], 'Step 2: Adjust/Draw Green Polygon to KEEP edges.'});

                    % 4. The user adjusts or draws the polygon; it is written into roilist
                    valid_selection = false;
                    while ~valid_selection && isvalid(hF)
                        try
                            if ~isempty(existing_pos)
                                poly = drawpolygon('Color', 'g', 'LineWidth', 2, 'Position', existing_pos);
                                existing_pos = []; % a redraw after a delete starts fresh
                            else
                                poly = drawpolygon('Color', 'g', 'LineWidth', 2);
                            end
                            wait(poly);

                            if isvalid(poly)
                                user_mask = createMask(poly);
                                vtx = poly.Position;
                                nowt = datetime('now');
                                if ~isempty(roi_idx)
                                    roilist.ROIs(roi_idx).Vertices = vtx;
                                    roilist.ROIs(roi_idx).Mask = user_mask;
                                    roilist.ROIs(roi_idx).Modified = nowt;
                                    roilist.ROIs(roi_idx).ImageSize = size(user_mask);
                                else
                                    newroi.Label     = contour_name;
                                    newroi.Mode      = 'polygon';
                                    newroi.Vertices  = vtx;
                                    newroi.Mask      = user_mask;
                                    newroi.ImageSize = size(user_mask);
                                    newroi.RefSlice  = 1;
                                    newroi.ROISlice  = [];
                                    newroi.Created   = nowt;
                                    newroi.Modified  = nowt;
                                    roilist.ROIs(end+1) = newroi;
                                end
                                disp(['Saved selection for ' disp_msg]);
                                valid_selection = true;
                            elseif isvalid(hF)
                                disp('Polygon deleted. Returning to drawing mode...');
                                title({['Step 1: Detected Edges for ' disp_msg], 'Polygon deleted. Please Draw New Polygon.'});
                            else
                                disp('Figure closed. Skipping...');
                            end
                        catch ME
                            if ~isvalid(hF)
                                disp('Figure closed during interaction.');
                                break;
                            end
                            disp(['Selection error: ' ME.message]);
                            break;
                        end
                    end

                    if isvalid(hF)
                        close(hF);
                    end
                end
            end
        end

        function transform(obj, roilist)
        % The six contours -> image-frame polar profiles, centred on the constricted BV
        % centroid; the edges are the ones drawcontour showed, kept by polygon. The four
        % constricted/dilated labels are required; an intermediate without one gives NaN
            state_of   = ["constricted", "constricted", "dilated", "dilated", "intermediate", "intermediate"];
            channel_of = ["bv", "pvs", "bv", "pvs", "bv", "pvs"];
            contour_labels = state_of + upper(channel_of) + "_contour";
            if isempty(roilist) || ~all(ismember(contour_labels(1:4), roilist.list()))
                warning('roilist is missing one of the 4 %s labels. Skipping polar transform.', 'contour');
                return;
            end

            % 1. Centre from the constricted BV mask, PAX angle for display
            try
                mask_cbv = roilist.getmask('constrictedBV_contour');
                if ~isempty(mask_cbv)
                    img_size = size(mask_cbv);
                    props = regionprops(mask_cbv, 'Centroid');
                    if ~isempty(props)
                        center = props(1).Centroid; % [x, y]
                    else
                        center = size(mask_cbv, [2, 1]) / 2;
                    end
                else
                    warning('Constricted BV mask found but empty. Using cluster data size for center.');
                    img_size = size(obj.clust_med_bv, [1, 2]);
                    center = [img_size(2), img_size(1)] / 2;
                end
            catch
                warning('Constricted BV mask not found in roilist. Cannot determine center.');
                return;
            end

            try
                pax_vertices = roilist.getvertices('pax');
                if ~isempty(pax_vertices) && size(pax_vertices, 1) >= 2
                    pax_vec = pax_vertices(2,:) - pax_vertices(1,:);
                    pax_deg = atan2d(pax_vec(2), pax_vec(1));
                else
                    warning('PAX ROI not found or invalid in roilist. Using 0 degrees.');
                    pax_deg = 0;
                end
            catch
                warning('Could not retrieve PAX vertices. Using 0 degrees.');
                pax_deg = 0;
            end

            % 2. Image-frame angle and radius maps; PAX is a display marker only
            frame = polar_frame([img_size(1), img_size(2)], center, 'n_sectors', 24);

            % 3. Per contour: edges inside the polygon -> farthest radius per 1-degree bin
            %    order matches the colours in analysis_polar_makefig: CBV CPVS DBV DPVS IBV IPVS
            theta_centers = deg2rad((0:359) + 0.5);

            profiles = struct('roi_label', {}, 'binned_r', {});
            for i = 1:numel(contour_labels)
                roi_label = char(contour_labels(i));
                try
                    roi_mask = roilist.getmask(roi_label);
                catch
                    roi_mask = [];
                end

                img = obj.stateimage(channel_of(i), state_of(i));
                try
                    edges = edge(imgaussfilt(double(img), 3), 'log', 0.005);
                catch
                    edges = false(size(img));
                end

                % Missing mask or no edge pixels -> all-NaN profile, so the array stays 6 long
                if isempty(roi_mask)
                    profiles(i) = struct('roi_label', roi_label, 'binned_r', nan(1, 360));
                    continue;
                end
                [rows, cols] = find(edges & roi_mask);
                if isempty(rows)
                    profiles(i) = struct('roi_label', roi_label, 'binned_r', nan(1, 360));
                    continue;
                end

                tr = roivtx2polar([cols, rows], frame.angle_map, frame.radius_map); % [rad; deg; radius]
                binned_r = bin_distal(tr(2, :), tr(3, :));
                profiles(i) = struct('roi_label', roi_label, 'binned_r', binned_r);
            end

            obj.polar_theta = theta_centers;
            obj.polar_profiles = profiles;
            obj.polar_binstart_deg = frame.bin_start_deg;
            obj.pax_angle = pax_deg;
        end

        function save2disk(obj, save_dir)
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
            polarcluster = obj;
            save(fullfile(save_dir, 'polarcluster.mat'), 'polarcluster');
        end
    end

    methods (Static)
        function obj = from_struct(s)
        % A polarcluster.mat saved as a struct; the intermediate is chosen here if it is missing
            obj = analysis_clusterpolar();
            obj.param.cluster_num = s.cluster_num;
            stored = ["kgph_lumen_columnidx", "clusterboundary", "clust_med_bv", "clust_med_csf", ...
                "cluster_summary", "filtered_clusteridx", "constricted_cluster_idx", ...
                "dilated_cluster_idx", "intermediate_cluster_idx", "polar_theta", ...
                "polar_profiles", "polar_binstart_deg", "pax_angle"];
            for f = stored
                if isfield(s, f)
                    obj.(f) = s.(f);
                end
            end
            if isempty(obj.intermediate_cluster_idx)
                obj.chooseintermediate();
            end
        end
    end

    methods (Access = private)
        function summarize(obj, upper_bv, lower_bv)
        % Per cluster: mean and std of the raw BV diameter, frame count, mean shift of each wall
            bv_diameter = lower_bv - upper_bv;
            low_change = lower_bv - median(lower_bv);
            up_change = median(upper_bv) - upper_bv;
            sorted_diameter = bv_diameter(obj.kgph_lumen_columnidx);
            sorted_low = low_change(obj.kgph_lumen_columnidx);
            sorted_up = up_change(obj.kgph_lumen_columnidx);

            obj.cluster_summary = struct('id', {}, 'mean_diameter', {}, 'std_diameter', {}, ...
                'duration', {}, 'mean_lowchange', {}, 'mean_upchange', {});
            for cluster_idx = 1:size(obj.clusterboundary, 1)
                span = obj.clusterboundary(cluster_idx, 1):obj.clusterboundary(cluster_idx, 2);
                obj.cluster_summary(cluster_idx).id = cluster_idx;
                obj.cluster_summary(cluster_idx).mean_diameter = mean(sorted_diameter(span));
                obj.cluster_summary(cluster_idx).std_diameter = std(sorted_diameter(span));
                obj.cluster_summary(cluster_idx).duration = numel(span);
                obj.cluster_summary(cluster_idx).mean_lowchange = mean(sorted_low(span));
                obj.cluster_summary(cluster_idx).mean_upchange = mean(sorted_up(span));
            end
        end

        function chooseintermediate(obj)
        % Nearest the median filtered diameter, among unpicked clusters at least as long as the
        % shorter pick; unweighted median, see CLAUDE_LOG.md "two states to three"
            cluster_diameter = [obj.cluster_summary.mean_diameter];
            cluster_duration = [obj.cluster_summary.duration];
            picked_idx = [obj.constricted_cluster_idx, obj.dilated_cluster_idx];
            frame_floor = min(cluster_duration(picked_idx));
            target_diameter = median(cluster_diameter(obj.filtered_clusteridx));
            is_long_enough = cluster_duration(obj.filtered_clusteridx) >= frame_floor;
            is_unpicked = ~ismember(obj.filtered_clusteridx, picked_idx);
            keep_candidate = is_long_enough & is_unpicked;
            candidate_idx = obj.filtered_clusteridx(keep_candidate);
            [~, nearest] = min(abs(cluster_diameter(candidate_idx) - target_diameter));
            obj.intermediate_cluster_idx = candidate_idx(nearest);
            fprintf('intermediate cluster : %d of %d candidates, %d frames, %.2f px against a median of %.2f\n', ...
                obj.intermediate_cluster_idx, numel(candidate_idx), ...
                cluster_duration(obj.intermediate_cluster_idx), ...
                cluster_diameter(obj.intermediate_cluster_idx), target_diameter);
        end
    end
end


function binned_r = bin_distal(angles_deg, radii)
% Farthest radius per 1-degree bin (0..359); an empty bin stays NaN so the line breaks there
deg_bins = 0:360;
binned_r = nan(1, 360);
for k = 1:360
    idx = angles_deg >= deg_bins(k) & angles_deg < deg_bins(k+1);
    if any(idx)
        binned_r(k) = max(radii(idx));
    end
end
end
