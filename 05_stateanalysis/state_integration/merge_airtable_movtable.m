function merged_table = merge_airtable_movtable(air_table, move_table, min_gap_sec)
% MERGE_AIRTABLE_MOVTABLE  Combine whisker-stim and movement time tables.
% Inputs
% air_table: Nx2 whisker-stim [start, end] (seconds)
% move_table: Mx2  [start, end]  ball-movement intervals (seconds)

%   Outputs
%     merged_table    : Kx2  merged [start, end] rows (sorted)
%     pure_move_table : Lx2  movement rows with no stim overlap and gap>=min_gap_sec

if nargin < 3 || isempty(min_gap_sec)
    min_gap_sec = 10;
end

n_air  = size(air_table,  1);
n_move = size(move_table, 1);

%%  1. Purify air puff without contaminated start by voluntary locomotion
% A "good" air puff is one where the mouse was NOT moving in the 10 seconds prior.
good_start = true(n_air, 1);

for a = 1:n_air
    ai_start = air_table(a, 1);
    for m = 1:n_move
        mv_start = move_table(m, 1);
        mv_end   = move_table(m, 2);
        
        % If a movement ends within 10 seconds before the air puff starts,
        % or if the mouse is already moving when the puff starts:
        if mv_end > ai_start - min_gap_sec && mv_start < ai_start
            good_start(a) = false;
            break; % Contaminated, skip checking other movements
        end
    end
end

filtered_air = air_table(good_start, :);
n_filteredair = size(filtered_air, 1);


%%  1. Purify air puff starting point
is_overlapping    = false(n_move, 1);   % overlaps an air interval
matched_air       = zeros(n_move, 1);   % which air row it overlaps (first)

for a = 1:n_filteredair
    ai_start = filtered_air(a, 1);
    ai_end   = filtered_air(a, 2);
    for m = 1:n_move
        mv_start = move_table(m, 1);
        mv_end   = move_table(m, 2);
        
        % The movement must start DURING the air puff (ai_start <= mv_start <= ai_end)
        % and it must continue to the end (or close to it)
        if mv_start >= ai_start && mv_start <= ai_end && mv_end + min_gap_sec >= ai_end
            is_overlapping(m) = true;
            matched_air(m)    = a;
        end
    end
end

% Step 2 – Build merged rows  (air start + max end)
merged_starts = filtered_air(:, 1);       % always use air start
merged_ends   = filtered_air(:, 2);       % initially set to air end
for m = 1:n_move
    if is_overlapping(m)
        a = matched_air(m);
        merged_ends(a) = max(merged_ends(a), move_table(m, 2));
    end
end
merged_table = [merged_starts, merged_ends];

end
