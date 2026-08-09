function out_path = piv_blink_gif(ch1, ch2, framesets, labels, out_path, opt)
%PIV_BLINK_GIF  Two averaged frames alternating, for eyeballing what moved.
%   The QC in front of PIV: if the blink shows no wall motion, no correlation
%   setting is going to find any. Each blink frame is the MEAN of a frame set, so
%   the two panels differ only in when they were taken.
%
% IN   ch1, ch2    H x W x T numeric  the two channels, whole recording
%                  ch2 = [] draws ch1 alone in grey
%      framesets   1 x 2 cell   frame indices averaged into panel A and panel B
%      labels      1 x 2 string burned into the top-left of each panel;
%                  "" on a panel draws no text there
%      out_path    char         .gif to write; the folder must already exist
%      cropmask    H x W bool   [] = whole frame, else crop to its bounding box
%      crop_pad    px around the bounding box, default 40
%      pixel2um    float um/px for the scale bar, default [] = no bar
%      scalebar_um float, default 20
%      delay_s     float s per panel, default 0.6
%      title_text  char         one line under the panels, default ''
% OUT  out_path
%
%   err   the two panels MUST share one contrast mapping. plot_make_rgb stretches
%         each image on its own percentiles, so a blink built from it shows the
%         brightness change between the two moments as if it were motion. The
%         limits here are taken once, over BOTH frame sets together
%   why   mean and not median : the sets are 10-50 frames of the same stable
%         state, so the mean is the higher-SNR estimator and there is no
%         transient for it to smear
%   caution  a GIF is 8-bit indexed. rgb2ind with 256 colours is enough for a
%            two-channel image but would band a smooth grey ramp

arguments
    ch1  {mustBeNumeric, mustBeNonempty}
    ch2  {mustBeNumeric} = []
    framesets   (1,2) cell = {}
    labels      (1,2) string = ["A", "B"]
    out_path    (1,:) char = ''
    opt.cropmask    logical = []
    opt.crop_pad    (1,1) double = 40
    opt.pixel2um    = []
    opt.scalebar_um (1,1) double = 20
    opt.delay_s     (1,1) double = 0.6
    opt.title_text  (1,:) char = ''
end

% 1. Average each frame set
panel = cell(1,2);
for p = 1:2
    idx = framesets{p};
    panel{p}.a = mean(double(ch1(:,:,idx)), 3);
    if ~isempty(ch2); panel{p}.b = mean(double(ch2(:,:,idx)), 3); else; panel{p}.b = []; end
end

% 2. ONE contrast mapping for both panels -- see the err note above
lim1 = prctile([panel{1}.a(:); panel{2}.a(:)], [10 99]);
if ~isempty(ch2)
    lim2 = prctile([panel{1}.b(:); panel{2}.b(:)], [10 99]);
end

% 3. Crop window from the mask's bounding box
if isempty(opt.cropmask)
    rows = 1:size(ch1,1);  cols = 1:size(ch1,2);
else
    [rr, cc] = find(opt.cropmask);
    rows = max(1, min(rr)-opt.crop_pad) : min(size(ch1,1), max(rr)+opt.crop_pad);
    cols = max(1, min(cc)-opt.crop_pad) : min(size(ch1,2), max(cc)+opt.crop_pad);
end

% 4. Render each panel through a figure so the bar is burned in
% note  an empty label or title draws nothing, and with no title the bottom strip
%       is not reserved at all -- the frame is then pure image plus scale bar
strip = 40 * ~isempty(opt.title_text);
fig_h = numel(rows)*2 + strip;
fig = figure('Visible', 'off', 'Color', 'k', 'Units', 'pixels', ...
    'Position', [100 100 numel(cols)*2 fig_h]);
frames = cell(1,2);
for p = 1:2
    clf(fig);
    ax = axes(fig, 'Position', [0 strip/fig_h 1 numel(rows)*2/fig_h]);
    a = mat2gray(panel{p}.a(rows, cols), double(lim1));
    if isempty(ch2)
        img = repmat(a, 1, 1, 3);
    else
        b = mat2gray(panel{p}.b(rows, cols), double(lim2));
        img = cat(3, a, b, a);          % ch1 magenta, ch2 cyan-green, as plot_make_rgb
    end
    image(ax, min(img, 1)); axis(ax, 'image', 'off');
    if strlength(labels(p)) > 0
        text(ax, 8, 18, labels(p), 'Color', 'w', 'FontSize', 15, ...
            'FontWeight', 'bold', 'FontName', 'Arial', 'Interpreter', 'none');
    end
    if ~isempty(opt.pixel2um)
        bar_px = opt.scalebar_um / opt.pixel2um;
        y = numel(rows) - 14;
        line(ax, [numel(cols)-14-bar_px, numel(cols)-14], [y y], ...
            'Color', 'w', 'LineWidth', 3);
    end
    if ~isempty(opt.title_text)
        annotation(fig, 'textbox', [0 0 1 strip/fig_h], ...
            'String', opt.title_text, 'Color', 'w', 'EdgeColor', 'none', ...
            'FontSize', 11, 'FontName', 'Arial', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    drawnow;
    frames{p} = frame2im(getframe(fig));
end
close(fig);

% 5. Write the two-panel loop
for p = 1:2
    [ind, cmap] = rgb2ind(frames{p}, 256, 'nodither');
    if p == 1
        imwrite(ind, cmap, out_path, 'gif', 'LoopCount', Inf, 'DelayTime', opt.delay_s);
    else
        imwrite(ind, cmap, out_path, 'gif', 'WriteMode', 'append', 'DelayTime', opt.delay_s);
    end
end
fprintf('blink gif  %s  (%d x %d px, %d + %d frames averaged)\n', out_path, ...
    size(frames{1},2), size(frames{1},1), numel(framesets{1}), numel(framesets{2}));
end
