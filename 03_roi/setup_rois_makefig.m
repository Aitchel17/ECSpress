function setup_rois_makefig(roilist, twophoton_processed, pax_fwhm, save_dir)
%SETUP_ROIS_MAKEFIG  Verification SVGs for the manual ROIs from setup_rois.
%   Each ROI is drawn over its OWN reference frame: addimgchannel rebuilds
%   ROISlice from RefSlice +/- 5, so a dilated ROI sits on the dilated frame.
%   Vessel ROIs go on the lumen channel and PVS ROIs on the PVS channel, both
%   read from the FWHM record because recordings swap ch1/ch2.
%
% IN   roilist              roi_handle carrying the manual_* ROIs
%      twophoton_processed  struct with ch1, ch2 (H x W x T)
%      pax_fwhm             line_fwhm; needs param.channel_lumen / channel_pvs
%      save_dir             where the SVGs are written
% OUT  none

arguments
    roilist
    twophoton_processed struct
    pax_fwhm
    save_dir (1,:) char
end

% 0. Setup
% 0.1 Channel slot showrois indexes in ROISlice: 1 = ch1, 2 = ch2, 3 = RGB of both
need = {'channel_lumen', 'channel_pvs'};
gone = need(~isfield(pax_fwhm.param, need));
if ~isempty(gone)
    error('setup_rois_makefig:channelMap', ...
        'pax_fwhm.param is missing %s.', strjoin(gone, ' and '));
end
bv_ch  = ch_index(pax_fwhm.param.channel_lumen);
pvs_ch = ch_index(pax_fwhm.param.channel_pvs);

% 0.2 Give every ROI its own background slice, taken around its RefSlice.
%     addimgchannel only acts on 4-D input, hence the cat along dim 4
ch1 = twophoton_processed.ch1;
ch2 = twophoton_processed.ch2;
if ndims(ch1) ~= 3 || ndims(ch2) ~= 3
    error('setup_rois_makefig:stackShape', 'ch1 and ch2 must be H x W x T.');
end
roilist.addimgchannel(cat(4, ch1, ch2), roilist.list());

% 1. Vessel boundaries on the lumen channel
make_one('roi_dilated_bv',      bv_ch,  ["manual_dilated_bv", "manual_constricted_bv"], ["-r", "-g"]);
make_one('roi_constricted_bv',  bv_ch,  ["manual_constricted_bv", "manual_dilated_bv"], ["-g", "-r"]);

% 2. PVS boundaries on the PVS channel
make_one('roi_dilated_pvs',     pvs_ch, ["manual_dilated_pvs", "manual_constricted_pvs"], ["-r", "-g"]);
make_one('roi_constricted_pvs', pvs_ch, ["manual_constricted_pvs", "manual_dilated_pvs"], ["-g", "-r"]);

% 3. Dip lines on the lumen channel
make_one('roi_dipinout',        bv_ch,  ["manual_dipin", "manual_dipout"], ["-g", "-r"]);

% 4. Axis, dip lines and vessel boundaries on the two-channel composite
make_one('roi_pax', 3, ...
    ["pax", "manual_dipin", "manual_dipout", "manual_constricted_bv", "manual_dilated_bv"], ...
    ["-c", "-y", "-y", "y", "y"]);

    function make_one(name, channel, labels, colors)
    %MAKE_ONE  One verification figure; a missing ROI skips its figure, not the run.
        try
            f = make_fig(name);
            f.update_figsize([8 6]);
            f.reset_axis();
            f.showrois(roilist, channel, labels, colors);
            f.save2svg(save_dir);
        catch ME
            warning(ME.identifier, 'Failed to generate %s: %s', name, ME.message);
        end
    end
end

% ---------------------------------------------------------------------------
function k = ch_index(name)
%CH_INDEX  'ch1' / 'ch2' -> the channel slot showrois indexes in ROISlice.
    k = str2double(erase(char(name), 'ch'));
    if ~ismember(k, [1 2])
        error('setup_rois_makefig:badChannel', ...
            'channel must be ch1 or ch2, got "%s".', char(name));
    end
end
