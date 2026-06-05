function transition_table = get_transition(transition_window,pre_timetable,post_timetable, tolerance)
if nargin < 4
    tolerance = 0;
end

% 1. Filter bouts by duration
if isempty(pre_timetable) || isempty(post_timetable)
    transition_table = [];
    return
end

pre_valid = pre_timetable( (pre_timetable(:,2) - pre_timetable(:,1)) >= transition_window, :);
post_valid = post_timetable( (post_timetable(:,2) - post_timetable(:,1)) >= transition_window, :);

%% Find Transitions
% Logic: Transition occurs if Ref Bout End Time == Post Bout Start Time
if isempty(pre_valid) || isempty(post_valid)
    trans_sec = [];
else
    if tolerance == 0
        translogic = ismember(pre_valid(:,2), post_valid(:,1));
        trans_sec  = pre_valid(translogic,2);
    else
        trans_sec = [];
        for i = 1:size(pre_valid, 1)
            % Check if any post event starts right after this pre event ends
            % (Allowing for overlap or gap up to 'tolerance' seconds)
            gaps = post_valid(:, 1) - pre_valid(i, 2);
            if any(abs(gaps) <= tolerance)
                trans_sec(end+1, 1) = pre_valid(i, 2); %#ok<AGROW>
            end
        end
    end
end

%%
transition_table = NaN([length(trans_sec),2]);
transition_table(:,1) = trans_sec - transition_window;
transition_table(:,2) = trans_sec + transition_window;

end


