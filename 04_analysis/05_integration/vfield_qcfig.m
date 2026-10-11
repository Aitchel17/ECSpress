function [quiver_fig, blink_fig] = vfield_qcfig(from_img, to_img, xyuv)
%VFIELD_QCFIG  One displacement field's QC pair : the field over the FROM image, and the FROM/TO blink.
%   Caller: analysis_main_piv
%   Both windows are made invisible and untitled; the caller titles them, copies their axes onto a
%   board and closes them.
%
% IN   from_img  H x W float        the FROM frames' mean : the quiver's background, the blink's first image
%      to_img    H x W float        the TO frames' mean, the blink's second image
%      xyuv      ny x nx x 4 float  the gated field, px; NaN = no vector
% OUT  quiver_fig  showpiv    arrows x arrow_scale, every other grid node drawn
%      blink_fig   make_fig   one image tagged blink_image, its UserData {from_rgb, to_rgb}; plot_blink
%                             flips it
    arrow_scale = 10;           % arrow length, times the px displacement
    contrast_pct = [10 99];     % one grey mapping for both images : percentiles of their pooled pixels

    % the field, every other grid row and column, stamped on the image grid
    keep = ~isnan(xyuv(:, :, 3)) & ~isnan(xyuv(:, :, 4));
    keep(2:2:end, :) = false;
    keep(:, 2:2:end) = false;
    uv_map = piv_stamp(xyuv, size(from_img), keep);
    quiver_fig = showpiv('vfield_qcfig quiver', arrow_scale);
    quiver_fig.fig.Visible = 'off';
    quiver_fig.plot_background(from_img);
    quiver_fig.plot_quiver(uv_map);

    % the blink : both images on one contrast, the second carried on the image for plot_blink
    grey_lim = prctile([from_img(:); to_img(:)], contrast_pct);
    from_rgb = repmat(mat2gray(from_img, grey_lim), [1 1 3]);
    to_rgb = repmat(mat2gray(to_img, grey_lim), [1 1 3]);
    blink_fig = make_fig('vfield_qcfig blink', 'normal');
    blink_fig.fig.Visible = 'off';
    image(blink_fig.ax, from_rgb, 'Tag', 'blink_image', 'UserData', {from_rgb, to_rgb});
    axis(blink_fig.ax, 'image');
end
