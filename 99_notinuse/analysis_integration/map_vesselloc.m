%% map_vesselloc.m
% Logic for mapping vessel locations across different imaging sessions
% based on relative positions of the objective lens.

clc; clear; close all;

%% 1. Generate Mock Data
% Simulate dirstruct (session_table)
dirstruct = generate_mock_data();

disp('Mock Data Generated:');
disp(head(dirstruct));

%% 2. Extract Data Per Date
dates = unique(dirstruct.date);
n_dates = numel(dates);
date_data = struct();

figure('Name', 'Raw Positions', 'Color', 'w');
hold on;
colors = lines(n_dates);
legend_entries = {};


for i = 1:n_dates
    current_date = dates{i};
    % Filter rows for this date
    rows = strcmp(dirstruct.date, current_date);
    sub_table = dirstruct(rows, :);
    
    % Extract XY positions
    xy = [cell2mat(sub_table.x), cell2mat(sub_table.y)];
    ids = sub_table.sess;
    
    % Store in struct
    date_data(i).Date = current_date;
    date_data(i).XY = xy;
    date_data(i).IDs = ids;
    
    % Plot raw data
    scatter(xy(:,1), xy(:,2), 50, colors(i,:), 'filled', 'DisplayName', current_date);
    text(xy(:,1), xy(:,2), ids, 'Color', colors(i,:), 'VerticalAlignment', 'bottom');
    legend_entries{end+1} = current_date;
end

title('Raw Objective Positions (Shifted)');
xlabel('X (\mum)'); ylabel('Y (\mum)');
grid on;
legend(legend_entries);
hold off;

%% 3. Match Across Dates using Relative Distances
% Algorithm:
% 1. Calculate pairwise distance matrix for specific date (Reference).
% 2. For Target Date, calculate pairwise distance matrix.
% 3. Find matching distances (within tolerance).
% 4. Identify corresponding nodes based on matching distance patterns.

ref_idx = 1; % Use first date as reference
ref_date = date_data(ref_idx);
ref_xy = ref_date.XY;
ref_ids = ref_date.IDs;

% Calculate Distance Matrix for Reference
% D_ref(i,j) is distance between point i and j
D_ref = squareform(pdist(ref_xy));

matches_table = table();

for i = 2:n_dates
    target_date = date_data(i);
    target_xy = target_date.XY;
    target_ids = target_date.IDs;
    
    % Calculate Distance Matrix for Target
    D_target = squareform(pdist(target_xy));
    
    fprintf('\nMatching %s -> %s\n', ref_date.Date, target_date.Date);
    
    tolerance = 10.0; % microns tolerance for distance matching
    
    % Store matches: [Res_Index, Target_Index]
    current_matches = [];
    
    for r = 1:size(ref_xy, 1)
        % Get distances from this ref point to all other ref points
        dists_r = D_ref(r, :);
        dists_r = dists_r(dists_r > 0); % exclude self-distance
        
        best_match_idx = -1;
        best_match_score = -1;
        
        for t = 1:size(target_xy, 1)
            % Get distances from this target point to all other target points
            dists_t = D_target(t, :);
            dists_t = dists_t(dists_t > 0);
            
            % Count how many distances match
            match_count = 0;
            for k = 1:length(dists_r)
                if any(abs(dists_t - dists_r(k)) < tolerance)
                    match_count = match_count + 1;
                end
            end
            
            % RELAXED THRESHOLD:
            % If we have very few points (e.g., 2), match_count can be 1.
            % So we accept match_count >= 1 if it's the best score.
            if match_count >= 1 && match_count > best_match_score
                best_match_score = match_count;
                best_match_idx = t;
            end
        end
        
        if best_match_idx > 0
            fprintf('  Match Found: %s <--> %s (Score: %d)\n', ...
                ref_ids{r}, target_ids{best_match_idx}, best_match_score);
            current_matches(end+1, :) = [r, best_match_idx];
        end
    end

    % Only keep unique matches (simple greedy approach might duplicate)
    if ~isempty(current_matches)
        [~, unique_idx] = unique(current_matches(:,2)); % unique target indices
        current_matches = current_matches(unique_idx, :);
    end
    
    % --- Visualization of Matches ---
    % Calculate translation based on first match to overlay
    if ~isempty(current_matches)
        p_ref = ref_xy(current_matches(1,1), :);
        p_tar = target_xy(current_matches(1,2), :);
        translation = p_ref - p_tar;
        
        aligned_target_xy = target_xy + translation;
        
        figure('Name', sprintf('Match: %s vs %s', ref_date.Date, target_date.Date));
        hold on;
        scatter(ref_xy(:,1), ref_xy(:,2), 100, 'b', 'filled', 'DisplayName', ref_date.Date);
        scatter(aligned_target_xy(:,1), aligned_target_xy(:,2), 100, 'r', 'filled', 'MarkerFaceAlpha', 0.5, 'DisplayName', [target_date.Date ' (Aligned)']);
        text(ref_xy(:,1), ref_xy(:,2), ref_ids, 'Color', 'b', 'VerticalAlignment', 'top');
        text(aligned_target_xy(:,1), aligned_target_xy(:,2), target_ids, 'Color', 'r', 'VerticalAlignment', 'bottom');
        
        % Draw lines for matches
        for m = 1:size(current_matches, 1)
            idx_r = current_matches(m,1);
            idx_t = current_matches(m,2);
            plot([ref_xy(idx_r,1), aligned_target_xy(idx_t,1)], ...
                 [ref_xy(idx_r,2), aligned_target_xy(idx_t,2)], 'k--');
        end
        
        legend;
        grid on;
        title(sprintf('Alignment based on %s-%s', ref_ids{current_matches(1,1)}, target_ids{current_matches(1,2)}));
        hold off;
    else
        fprintf('  No valid matches found to align %s.\n', target_date.Date);
    end
end


%% Helper Function: Mock Data Generator
function df = generate_mock_data()
    % Create base positions (Date 1) - GROUND TRUTH 5 VESSELS
    true_xy = [0, 0;
               0, 150; 
               200, 200; 
               -50, 300]; 
    
    n_xy = size(true_xy, 1);
    
    % Define dates
    date_list = {'1', '2', '3'};
    
    % Initialize cell arrays for columns
    col.date = {};
    col.sess = {};
    col.x = {};
    col.y = {};
    col.z = {};
    
    % --- Date 1: Random selection of 2-3 points ---
    n_session = randi([3, 4]); % Select 2 or 3
    idx_1 = sort(randperm(n_xy, n_session)); % from true_xy select # n_session of idx
    
    % Apply Identity Transform
    xy_1 = true_xy(idx_1, :);
    
    for i = 1:length(idx_1)
        col.date{end+1,1} = date_list{1};
        col.sess{end+1,1} = sprintf('sess_%02d', idx_1(i)); % Maintain ID for verification
        col.x{end+1,1} = xy_1(i,1);
        col.y{end+1,1} = xy_1(i,2);
        col.z{end+1,1} = 50;
    end
    
    % --- Date 2: Different Random Subset + Shift + Noise ---
    n_session = randi([3, 4]);
    idx_2 = sort(randperm(n_xy, n_session));
    
    shift = [500, 500];
    xy_2 = true_xy(idx_2, :) + shift + randn(n_session, 2)*2;
    
    for i = 1:length(idx_2)
        col.date{end+1,1} = date_list{2};
        col.sess{end+1,1} = sprintf('sess_%02d', idx_2(i));
        col.x{end+1,1} = xy_2(i,1);
        col.y{end+1,1} = xy_2(i,2);
        col.z{end+1,1} = 55;
    end

    % --- Date 3: Another Random Subset + Shift + Rotation ---
    n_session = randi([2, 3]);
    idx_3 = sort(randperm(n_xy, n_session));
    
    shift3 = [-200, 1000];
    
    xy_3 = true_xy(idx_3, :) + shift3 + randn(n_session, 2)*2;
    
    for i = 1:length(idx_3)
        col.date{end+1,1} = date_list{3};
        col.sess{end+1,1} = sprintf('sess_%02d', idx_3(i));
        col.x{end+1,1} = xy_3(i,1);
        col.y{end+1,1} = xy_3(i,2);
        col.z{end+1,1} = 45;
    end
    
    % Construct Table
    df = struct2table(col);
end
