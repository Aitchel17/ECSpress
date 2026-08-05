function frame_range = opticalflow_viewer(imgstack, flow, opt)
%OPTICALFLOW_VIEWER  Interactive frame-by-frame quiver viewer.
%   Draws one flow pair at a time over its background image. Arrow buttons or
%   arrow keys step, [Mark Start] / [Mark End] set the range. Blocks until the
%   figure is closed; every index here is a flow-pair index, not an image frame.
%
% IN   imgstack     HxWxN preprocessed image stack
%      flow         HxWxNpairs x2 field from opticalflow_wlet
%      block_size   spatial averaging block size (px), default 10
%      scale        quiver arrow scale factor, default 1
%      uniform      normalize arrows to unit length, default false
%      color        RGB quiver color, default [0 1 0]
%      linewidth    arrow line width, default 0.5
% OUT  frame_range  [start end] pair indices marked by the user, default
%                   [1 n_pairs] if never marked

arguments
    imgstack         {mustBeNumeric, mustBeNonempty}
    flow             {mustBeNumeric, mustBeNonempty}
    opt.block_size   double  = 10
    opt.scale        double  = 1
    opt.uniform      logical = false
    opt.color        double  = [0 1 0]
    opt.linewidth    double  = 0.5
end

% 0. Setup
% 0.1 Pair count
n_pairs = size(flow, 3);
if n_pairs < 1
    error('opticalflow_viewer: flow has no pairs.');
end
% 0.2 Range state, shared with the callbacks
sel_start = 1;
sel_end   = n_pairs;

% 1. Build the figure
% 1.1 Figure window
fig = figure('Name', 'Optical Flow Viewer', ...
    'NumberTitle', 'off', ...
    'Color', [0.15 0.15 0.15], ...
    'KeyPressFcn', @on_key, ...
    'CloseRequestFcn', @on_close);

% 1.2 Image axes
ax = axes('Parent', fig, 'Position', [0.02 0.16 0.96 0.82]);

% 1.3 Pair slider
uicontrol('Style', 'slider', ...
    'Units', 'normalized', 'Position', [0.08 0.07 0.84 0.04], ...
    'Min', 1, 'Max', n_pairs, 'Value', 1, ...
    'SliderStep', [1/(n_pairs-1+eps)  5/(n_pairs-1+eps)], ...
    'Callback', @on_slider, ...
    'Tag', 'frame_slider');

% 1.4 Prev / next buttons
uicontrol('Style', 'pushbutton', 'String', char(9664), ...
    'Units', 'normalized', 'Position', [0.01 0.07 0.06 0.04], ...
    'Callback', @(~,~) step(-1));
uicontrol('Style', 'pushbutton', 'String', char(9654), ...
    'Units', 'normalized', 'Position', [0.93 0.07 0.06 0.04], ...
    'Callback', @(~,~) step(+1));

% 1.5 Mark start / mark end buttons
uicontrol('Style', 'pushbutton', 'String', 'Mark Start', ...
    'Units', 'normalized', 'Position', [0.01 0.01 0.15 0.05], ...
    'BackgroundColor', [0.2 0.5 0.2], 'ForegroundColor', 'w', ...
    'Callback', @on_mark_start);
uicontrol('Style', 'pushbutton', 'String', 'Mark End', ...
    'Units', 'normalized', 'Position', [0.17 0.01 0.15 0.05], ...
    'BackgroundColor', [0.5 0.2 0.2], 'ForegroundColor', 'w', ...
    'Callback', @on_mark_end);

% 1.6 Marked-range label
lbl_range = uicontrol('Style', 'text', ...
    'Units', 'normalized', 'Position', [0.33 0.01 0.35 0.05], ...
    'BackgroundColor', [0.15 0.15 0.15], 'ForegroundColor', [0.8 0.8 0.8], ...
    'FontSize', 9, 'Tag', 'range_label');

% 1.7 Pair counter label
lbl_frame = uicontrol('Style', 'text', ...
    'Units', 'normalized', 'Position', [0.30 0.12 0.40 0.03], ...
    'BackgroundColor', [0.15 0.15 0.15], 'ForegroundColor', 'w', ...
    'FontSize', 10, 'Tag', 'frame_label');

% 2. Show the first pair
update_range_label();
draw_frame(1);

% 3. Block until the figure closes, then hand back the marked range
waitfor(fig);
frame_range = [sel_start, sel_end];

% 4. Callbacks
    function on_slider(src, ~)
        %ON_SLIDER  Draw the pair the slider points at.
        draw_frame(round(src.Value));
    end

    function on_key(~, evt)
        %ON_KEY  Step one pair on an arrow key.
        switch evt.Key
            case {'rightarrow', 'downarrow'};  step(+1);
            case {'leftarrow',  'uparrow'};    step(-1);
        end
    end

    function step(d)
        %STEP  Move d pairs from the slider position and redraw.
        %   The target index is clamped to [1 n_pairs].
        sl = findobj(fig, 'Tag', 'frame_slider');
        k  = max(1, min(n_pairs, round(sl.Value) + d));
        sl.Value = k;
        draw_frame(k);
    end

    function on_mark_start(~, ~)
        %ON_MARK_START  Mark the current pair as the range start.
        %   Drags the end along if it would sit before the new start.
        sel_start = round(findobj(fig,'Tag','frame_slider').Value);
        if sel_start > sel_end; sel_end = sel_start; end
        update_range_label();
    end

    function on_mark_end(~, ~)
        %ON_MARK_END  Mark the current pair as the range end.
        %   Drags the start back if it would sit after the new end.
        sel_end = round(findobj(fig,'Tag','frame_slider').Value);
        if sel_end < sel_start; sel_start = sel_end; end
        update_range_label();
    end

    function on_close(~, ~)
        %ON_CLOSE  Destroy the figure.
        delete(fig);   % releases waitfor
    end

    function update_range_label()
        %UPDATE_RANGE_LABEL  Refresh the marked-range text.
        lbl_range.String = sprintf('Range: %d  →  %d', sel_start, sel_end);
    end

% 5. Frame drawing
    function draw_frame(k)
        %DRAW_FRAME  Draw pair k: background image, quiver overlay, counter.
        % 5.1 Background image
        bg = imgstack(:,:,k);
        axes(ax); %#ok<LAXES>
        imagesc(ax, bg); colormap(ax, gray); axis(ax, 'image');
        ax.XColor = 'w'; ax.YColor = 'w';
        hold(ax, 'on');

        % 5.2 Flow components on the pixel grid
        flow_k = squeeze(flow(:,:,k,:));
        u = flow_k(:,:,1);
        v = flow_k(:,:,2);
        [H, W] = size(u);
        [x, y] = meshgrid(1:W, 1:H);

        % 5.3 Drop magnitude, keep direction
        if opt.uniform
            mag = sqrt(u.^2 + v.^2); mag(mag==0) = 1;
            u = u./mag; v = v./mag;
        end

        % 5.4 Block-average onto a coarser arrow grid
        s = round(opt.block_size);
        if s > 1
            % 5.5 Crop to a whole number of s x s blocks
            Hr = floor(H/s)*s; Wr = floor(W/s)*s;
            u = u(1:Hr,1:Wr);  v = v(1:Hr,1:Wr);
            % 5.6 Average down rows, then across columns
            u = squeeze(mean(reshape(u,  s, Hr/s, Wr), 1));
            v = squeeze(mean(reshape(v,  s, Hr/s, Wr), 1));
            u = squeeze(mean(reshape(u', s, Wr/s, Hr/s), 1))';
            v = squeeze(mean(reshape(v', s, Wr/s, Hr/s), 1))';
            % 5.7 Anchor the arrows at block centres
            [x, y] = meshgrid((0.5:Wr/s)*s, (0.5:Hr/s)*s);
        end

        % 5.8 Quiver overlay
        quiver(ax, x, y, u*opt.scale, v*opt.scale, ...
            'Color', opt.color, 'AutoScale', 'off', 'LineWidth', opt.linewidth);
        hold(ax, 'off');

        % 5.9 Pair counter
        lbl_frame.String = sprintf('Pair %d / %d', k, n_pairs);
    end

end
