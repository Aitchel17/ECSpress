function setup_rois(roilist, twophoton_processed, pax_fwhm, primary_path, opt)
%SETUP_ROIS  Trace the BV and PVS boundaries on the FWHM diameter extremes.
%   A dilated or constricted boundary only exists at one frame each: where the
%   FWHM lumen trace peaks and where it bottoms. Each ROI is drawn on a short
%   window around its own frame rather than by scrolling the whole recording.
%   The diameter is stationary at an extreme, so the window averages cleanly.
%
% IN   roilist              roi_handle, modified in place
%      twophoton_processed  struct with ch1, ch2 (H x W x T)
%      pax_fwhm             line_fwhm; needs thickness.bv, param.fs and
%                           param.channel_lumen / channel_pvs
%      primary_path         roilist save location
%      halfwin              frames either side of an extreme, default 15
%      sg_win_s             smoothing window (s) for locating the extremes, default 3
% OUT  none; roilist is a handle, and is saved to primary_path

arguments
    roilist
    twophoton_processed struct
    pax_fwhm
    primary_path (1,:) char
    opt.halfwin  (1,1) double {mustBePositive} = 15
    opt.sg_win_s (1,1) double {mustBePositive} = 3
end

% 0. Setup
% 0.1 Channel map, read from the FWHM record: recordings swap ch1/ch2
need = {'channel_lumen', 'channel_pvs'};
gone = need(~isfield(pax_fwhm.param, need));
if ~isempty(gone)
    error('setup_rois:channelMap', ...
        ['pax_fwhm.param is missing %s. Re-run the FWHM cells of analysis_main, ' ...
         'which set them from bvch / csfch.'], strjoin(gone, ' and '));
end
bvimg  = twophoton_processed.(pax_fwhm.param.channel_lumen);
pvsimg = twophoton_processed.(pax_fwhm.param.channel_pvs);

% 0.2 Widest and narrowest frames, located on a smoothed lumen trace
dtrace  = double(pax_fwhm.thickness.bv(:)).';
sg_win  = 2*floor(opt.sg_win_s * pax_fwhm.param.fs / 2) + 1;
dsg     = sgolayfilt(dtrace, 3, sg_win);
[~, f_dilated]     = max(dsg);
[~, f_constricted] = min(dsg);

% 0.3 Frame window around each extreme. pre_groupaverage discards an incomplete
%     group of 10, so the window must stay well above 10 frames
nT    = size(bvimg, 3);
w_dil = max(1, f_dilated     - opt.halfwin) : min(nT, f_dilated     + opt.halfwin);
w_con = max(1, f_constricted - opt.halfwin) : min(nT, f_constricted + opt.halfwin);
report_frame('dilated',     f_dilated,     dtrace, pax_fwhm);
report_frame('constricted', f_constricted, dtrace, pax_fwhm);

% 1. PVS boundaries, drawn on the PVS channel
draw_step(roilist, pvsimg, bvimg, 'manual_dilated_pvs',     '',                       'polygon', w_dil, primary_path);
draw_step(roilist, pvsimg, bvimg, 'manual_constricted_pvs', 'manual_dilated_pvs',     'polygon', w_con, primary_path);

% 2. Vessel boundaries, drawn on the lumen channel
draw_step(roilist, bvimg, pvsimg, 'manual_dilated_bv',      'manual_dilated_pvs',     'polygon', w_dil, primary_path);
draw_step(roilist, bvimg, pvsimg, 'manual_constricted_bv',  'manual_constricted_pvs', 'polygon', w_con, primary_path);

% 3. Dip lines, on the dilated frame where the structures are widest
draw_step(roilist, bvimg,  pvsimg, 'manual_dipin',  '', 'line', w_dil, primary_path);
draw_step(roilist, pvsimg, bvimg,  'manual_dipout', '', 'line', w_dil, primary_path);
end

% ---------------------------------------------------------------------------
function report_frame(name, f, dtrace, pax_fwhm)
%REPORT_FRAME  Print which frame an extreme landed on.
    if isprop(pax_fwhm, 't_axis') && numel(pax_fwhm.t_axis) >= f
        t = pax_fwhm.t_axis(f);
    else
        t = f / pax_fwhm.param.fs;
    end
    fprintf('%-12s frame %5d   t = %7.1f s   bv = %.2f px\n', name, f, t, dtrace(f));
end

% ---------------------------------------------------------------------------
function draw_step(roilist, refimg, otherimg, label, source, mode, w, primary_path)
%DRAW_STEP  Create, copy-then-adjust, or modify one ROI on the frame window w.
%   refimg is the channel shown; otherimg is the co-registered reflect view.
%   Saves as soon as the ROI is finished: it is hand-drawn and cannot be redone,
%   and a later step failing must not cost the ones already traced.
    % 0. The two windowed stacks handed to the selector
    ref   = refimg(:, :, w);
    other = otherimg(:, :, w);

    % 1. Existing ROI: modify only on request
    if any(strcmp(roilist.list(), label))
        fprintf('\nROI "%s" already exists.\n', label);
        if ~strcmpi(input(sprintf('MODIFY "%s"? [y/N]: ', label), 's'), 'y')
            return
        end
        roilist.modifyroi(ref, label, other);
    else
        % 2. New ROI: seed from the source ROI when there is one
        fprintf('\nCreating NEW ROI: "%s"...\n', label);
        if ~isempty(source) && any(strcmp(roilist.list(), source))
            roilist.copyroi(source, label);
            fprintf('  -> copied from "%s", opening for adjustment\n', source);
            roilist.modifyroi(ref, label, other);
        else
            fprintf('  -> drawing from scratch\n');
            roilist.addormodifyroi(ref, label, mode, other);
        end
    end

    % 3. RefSlice comes back indexed against w; put it back on the full stack.
    %    findLabel is restricted to make_fig, so the search runs on ROIs directly
    i = find(strcmp({roilist.ROIs.Label}, label), 1);
    roilist.ROIs(i).RefSlice = w(1) - 1 + roilist.ROIs(i).RefSlice;

    % 4. Persist this ROI before moving to the next one
    roilist.save2disk(primary_path);
    fprintf('  -> saved (RefSlice %d)\n', roilist.ROIs(i).RefSlice);
end
