function frame_range = opticalflow_viewer(imgstack, flow, opt)
%OPTICALFLOW_VIEWER  Interactive frame-by-frame quiver viewer.
%   Displays the preprocessed image with optical flow overlay.
%   Use ◀/▶ buttons or ←/→ keys to navigate. Mark a start/end frame
%   with the [Mark Start] / [Mark End] buttons, then close the figure
%   to return the selected range.
%
%   USAGE
%     frame_range = opticalflow_viewer(imgstack, flow)
%     frame_range = opticalflow_viewer(imgstack, flow, block_size=15, scale=3)
%
%   INPUTS
%     imgstack   : double H x W x N array (preprocessed, background)
%     flow       : double H x W x N_pairs x 2 from opticalflow_wlet
%     block_size : spatial averaging block size (default: 10)
%     scale      : quiver arrow scale factor (default: 1)
%     uniform    : normalize to unit length (default: false)
%     color      : RGB quiver color (default: [0 1 0] green)
%     linewidth  : arrow line width (default: 0.5)
%
%   OUTPUT
%     frame_range : [start end] pair indices marked by the user
%                   (defaults to [1 n_pairs] if not explicitly marked)

arguments
    imgstack         {mustBeNumeric, mustBeNonempty}
    flow             {mustBeNumeric, mustBeNonempty}
    opt.block_size   double  = 10
    opt.scale        double  = 1
    opt.uniform      logical = false
    opt.color        double  = [0 1 0]
    opt.linewidth    double  = 0.5
end

n_pairs = size(flow, 3);
if n_pairs < 1
    error('opticalflow_viewer: flow has no pairs.');
end

% Shared range state
sel_start = 1;
sel_end   = n_pairs;

%% Build figure
fig = figure('Name', 'Optical Flow Viewer', ...
    'NumberTitle', 'off', ...
    'Color', [0.15 0.15 0.15], ...
    'KeyPressFcn', @on_key, ...
    'CloseRequestFcn', @on_close);

ax = axes('Parent', fig, 'Position', [0.02 0.16 0.96 0.82]);

% Slider
uicontrol('Style', 'slider', ...
    'Units', 'normalized', 'Position', [0.08 0.07 0.84 0.04], ...
    'Min', 1, 'Max', n_pairs, 'Value', 1, ...
    'SliderStep', [1/(n_pairs-1+eps)  5/(n_pairs-1+eps)], ...
    'Callback', @on_slider, ...
    'Tag', 'frame_slider');

% Prev / Next buttons
uicontrol('Style', 'pushbutton', 'String', char(9664), ...
    'Units', 'normalized', 'Position', [0.01 0.07 0.06 0.04], ...
    'Callback', @(~,~) step(-1));
uicontrol('Style', 'pushbutton', 'String', char(9654), ...
    'Units', 'normalized', 'Position', [0.93 0.07 0.06 0.04], ...
    'Callback', @(~,~) step(+1));

% Mark Start / Mark End buttons
uicontrol('Style', 'pushbutton', 'String', 'Mark Start', ...
    'Units', 'normalized', 'Position', [0.01 0.01 0.15 0.05], ...
    'BackgroundColor', [0.2 0.5 0.2], 'ForegroundColor', 'w', ...
    'Callback', @on_mark_start);
uicontrol('Style', 'pushbutton', 'String', 'Mark End', ...
    'Units', 'normalized', 'Position', [0.17 0.01 0.15 0.05], ...
    'BackgroundColor', [0.5 0.2 0.2], 'ForegroundColor', 'w', ...
    'Callback', @on_mark_end);

% Range label
lbl_range = uicontrol('Style', 'text', ...
    'Units', 'normalized', 'Position', [0.33 0.01 0.35 0.05], ...
    'BackgroundColor', [0.15 0.15 0.15], 'ForegroundColor', [0.8 0.8 0.8], ...
    'FontSize', 9, 'Tag', 'range_label');

% Frame label
lbl_frame = uicontrol('Style', 'text', ...
    'Units', 'normalized', 'Position', [0.30 0.12 0.40 0.03], ...
    'BackgroundColor', [0.15 0.15 0.15], 'ForegroundColor', 'w', ...
    'FontSize', 10, 'Tag', 'frame_label');

update_range_label();
draw_frame(1);

% Block until user closes the figure
waitfor(fig);
frame_range = [sel_start, sel_end];

%% ── Callbacks ──────────────────────────────────────────────────────────────
    function on_slider(src, ~)
        draw_frame(round(src.Value));
    end

    function on_key(~, evt)
        switch evt.Key
            case {'rightarrow', 'downarrow'};  step(+1);
            case {'leftarrow',  'uparrow'};    step(-1);
        end
    end

    function step(d)
        sl = findobj(fig, 'Tag', 'frame_slider');
        k  = max(1, min(n_pairs, round(sl.Value) + d));
        sl.Value = k;
        draw_frame(k);
    end

    function on_mark_start(~, ~)
        sel_start = round(findobj(fig,'Tag','frame_slider').Value);
        if sel_start > sel_end; sel_end = sel_start; end
        update_range_label();
    end

    function on_mark_end(~, ~)
        sel_end = round(findobj(fig,'Tag','frame_slider').Value);
        if sel_end < sel_start; sel_start = sel_end; end
        update_range_label();
    end

    function on_close(~, ~)
        delete(fig);   % releases waitfor
    end

    function update_range_label()
        lbl_range.String = sprintf('Range: %d  →  %d', sel_start, sel_end);
    end

%% ── Draw one frame ──────────────────────────────────────────────────────
    function draw_frame(k)
        bg = imgstack(:,:,k);
        axes(ax); %#ok<LAXES>
        imagesc(ax, bg); colormap(ax, gray); axis(ax, 'image');
        ax.XColor = 'w'; ax.YColor = 'w';
        hold(ax, 'on');

        flow_k = squeeze(flow(:,:,k,:));
        u = flow_k(:,:,1);
        v = flow_k(:,:,2);
        [H, W] = size(u);
        [x, y] = meshgrid(1:W, 1:H);

        if opt.uniform
            mag = sqrt(u.^2 + v.^2); mag(mag==0) = 1;
            u = u./mag; v = v./mag;
        end

        s = round(opt.block_size);
        if s > 1
            Hr = floor(H/s)*s; Wr = floor(W/s)*s;
            u = u(1:Hr,1:Wr);  v = v(1:Hr,1:Wr);
            u = squeeze(mean(reshape(u,  s, Hr/s, Wr), 1));
            v = squeeze(mean(reshape(v,  s, Hr/s, Wr), 1));
            u = squeeze(mean(reshape(u', s, Wr/s, Hr/s), 1))';
            v = squeeze(mean(reshape(v', s, Wr/s, Hr/s), 1))';
            [x, y] = meshgrid((0.5:Wr/s)*s, (0.5:Hr/s)*s);
        end

        quiver(ax, x, y, u*opt.scale, v*opt.scale, ...
            'Color', opt.color, 'AutoScale', 'off', 'LineWidth', opt.linewidth);
        hold(ax, 'off');

        lbl_frame.String = sprintf('Pair %d / %d', k, n_pairs);
    end

end
