function flipper = plot_blink(fig_obj, tag_name)
%PLOT_BLINK  Flip every image on a figure that carries tag_name between its two frames, until the figure closes.
%   Caller: a_analysis_piv
%
% IN   fig_obj      1x1 figure   the board
%      tag_name      1 x n char   the Tag of the images to flip; each one's UserData is {rgb_1, rgb_2}
% OUT  flipper  1x1 timer    running; the figure's DeleteFcn stops and deletes it
    period_s = 0.6;                    % time between flips (s)
    image_list = findobj(fig_obj, 'Tag', tag_name);
    flipper = timer('Period', period_s, 'ExecutionMode', 'fixedRate', 'TimerFcn', @flip_frame, ...
        'UserData', struct('image_list', image_list, 'phase', 1));
    fig_obj.DeleteFcn = @(~, ~) stop_blink(flipper);
    start(flipper);
end

function flip_frame(flipper, ~)
%FLIP_FRAME  One tick : every image to its other frame.
    tick = flipper.UserData;
    tick.phase = 3 - tick.phase;
    for k = 1:numel(tick.image_list)
        tick.image_list(k).CData = tick.image_list(k).UserData{tick.phase};
    end
    flipper.UserData = tick;
end

function stop_blink(flipper)
%STOP_BLINK  The figure closed : its timer goes with it.
    stop(flipper);
    delete(flipper);
end
