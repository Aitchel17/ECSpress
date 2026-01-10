function polarcluster = analysis_clusterpolar_contour(polarcluster, roilist)
% ANALYSIS_PAX_CLUSTER_CONTOUR Semi-automatic contouring of cluster median images.
%   Now integrates with roilist for persistence of the Selection Polygon.
%
%   1. Check roilist for existing 'contour_name' ROI.
%   2. Auto-detect edges.
%   3. Display edges + (Existing Polygon OR New Polygon).
%   4. User adjusts/draws polygon to KEEP edges.
%   5. Updates roilist with the Polygon (Vertices/Mask).
%   6. Updates polarcluster with the Filtered EdgeMask.
%
% Inputs:
%   polarcluster: Struct containing cluster images
%   roilist: roi_handle object
%
% Outputs:
%   polarcluster: Updated struct with 'manual_roi' containing filtered EdgeMask.

% Check if 4D data exists
if ~isfield(polarcluster, 'cluster_bvcsf_constdil')
    error('polarcluster must contain cluster_bvcsf_constdil (x,y,ch,state).');
end

% Initialize manual_roi if not exists
if ~isfield(polarcluster, 'manual_roi')
    polarcluster.manual_roi = struct();
end

state_names = {'constricted', 'dilated'};
ch_names    = {'bv', 'pvs'};

% Counter for linear indexing (1 to 4)
counter = 0;

% Iterate through States (1=Constricted, 2=Dilated) and Channels (1=BV, 2=PVS)
for state_idx = 1:2
    for ch_idx = 1:2
        counter = counter + 1;

        % Generate names and description
        contour_name = sprintf('%s%s_contour', state_names{state_idx}, upper(ch_names{ch_idx}));
        disp_msg     = sprintf('%s %s', state_names{state_idx}, ch_names{ch_idx});

        % Extract Image: (x, y, ch, state)
        img = polarcluster.cluster_bvcsf_constdil(:,:,ch_idx,state_idx);
        img_norm = mat2gray(img);

        % 1. Auto-detect edges
        disp(['Auto-detecting edges for: ' disp_msg '...']);
        try
            img_smooth = imgaussfilt(double(img), 3);
            edges = edge(img_smooth, 'log', 0.005);
        catch ME
            disp(['Error in auto-detection: ' ME.message]);
            edges = false(size(img));
        end

        % 2. Check ROILIST for existing polygon
        existing_pos = [];
        roi_idx = [];
        if ~isempty(roilist) && ~isempty(roilist.ROIs)
            roi_idx = find(strcmp({roilist.ROIs.Label}, contour_name), 1);
            if ~isempty(roi_idx)
                existing_pos = roilist.ROIs(roi_idx).Vertices;
                disp(['Found existing ROI in roilist: ' contour_name]);
            end
        end

        % 3. Open Figure for Selection
        hF = figure('Name', ['Select Edges: ' disp_msg], 'NumberTitle', 'off');
        imshow(img_norm, []);
        hold on;

        % Overlay detected edges in Red
        [Ay, Ax] = find(edges);
        if ~isempty(Ax)
            plot(Ax, Ay, 'r.', 'MarkerSize', 2);
        else
            text(10, 10, 'No edges detected', 'Color', 'y');
        end

        title({['Step 1: Detected Edges for ' disp_msg], 'Step 2: Adjust/Draw Green Polygon to KEEP edges.'});

        % 4. User Interaction
        valid_selection = false;

        while ~valid_selection && isvalid(hF)
            try
                if ~isempty(existing_pos)
                    poly = drawpolygon('Color', 'g', 'LineWidth', 2, 'Position', existing_pos);
                    existing_pos = []; % Reset so consecutive redraws start fresh if user deleted it
                else
                    poly = drawpolygon('Color', 'g', 'LineWidth', 2);
                end

                % Wait for completion
                custom_wait(poly);

                if isvalid(poly)
                    user_mask = createMask(poly);
                    vtx = poly.Position;

                    % 5. Filter Edges using Polygon
                    % filtered_edges = edges & user_mask; % Unused now

                    % 6. Update ROILIST (Persistence)
                    nowt = datetime('now');
                    if ~isempty(roi_idx)
                        % Modify existing
                        roilist.ROIs(roi_idx).Vertices = vtx;
                        roilist.ROIs(roi_idx).Mask = user_mask;
                        roilist.ROIs(roi_idx).Modified = nowt;
                        roilist.ROIs(roi_idx).ImageSize = size(user_mask);
                    else
                        % Create new
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
                else
                    if isvalid(hF)
                        disp('Polygon deleted. Returning to drawing mode...');
                        title({['Step 1: Detected Edges for ' disp_msg], 'Polygon deleted. Please Draw New Polygon.'});
                    else
                        disp('Figure closed. Skipping...');
                    end
                end

            catch ME
                if ~isvalid(hF)
                    disp('Figure closed during interaction.');
                    break;
                end
                disp(['Selection error: ' ME.message]);
                break; % Exit loop on error to prevent infinite loop
            end
        end



        if isvalid(hF)
            close(hF);
        end

    end
end

end

function custom_wait(roi)
% Helper to wait for ROI completion
wait(roi);
end
