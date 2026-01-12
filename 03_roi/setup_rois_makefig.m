function setup_rois_makefig(roilist, ch1, ch2, save_dir)
% SETUP_ROIS_MAKEFIG Generates verification figures for manual ROIs.
%   Re-implementation based on logic from legacy main_primary.m (lines 124-182).
%   Uses the make_fig class to generate and save SVG figures for validation.
%
%   Inputs:
%       roilist  - roi_handle object
%       ch1      - Ch1 image stack (for reference, though make_fig might use stored slices)
%       ch2      - Ch2 image stack
%       save_dir - Directory to save the figures

% Define the set of figures to generate, matching legacy logic
%% Update ROISlice with multi-channel data
% Concatenate ch1 and ch2 along the 4th dimension to create (H,W,T,C)
% roi_handle.addimgchannel requires 4D input to average over T (dim 3)
if ndims(ch1) == 3 && ndims(ch2) == 3
    rgb_stack = cat(4, ch1, ch2);
    labels = roilist.list();
    roilist.addimgchannel(rgb_stack, labels);
else
    warning('Input stacks are not 3D. Skipping addimgchannel update.');
end

% Define the set of figures to generate, matching legacy logic
%%
% here we need to make rgb stack by concatenating ch1 and ch2
rgb_stack = cat(3, ch1, ch2);
% This will be used to update roilist.ROIs.ROISlice which is required for showrois, without update ROIslice is just one channel and showrois channel will crush
labels = roilist.list;

roilist.addimgchannel(rgb_stack, labels)
%%


%% 5.1.1 BV_constricted + dilated and constricted vessel outline
try
    fig_dilatedbv = make_fig('roi_dilated_bv');
    fig_dilatedbv.update_figsize([8 6]);
    fig_dilatedbv.reset_axis();
    % Channel 1 background, red dilated, green constricted
    fig_dilatedbv.showrois(roilist, 1, ["manual_dilated_bv", "manual_constricted_bv"], ["-r", "-g"]);
    fig_dilatedbv.save2svg(save_dir);
catch ME
    warning(ME.identifier, 'Failed to generate roi_dilated_bv: %s', ME.message);
end

%% 5.1.2 BV_constricted + dilated and constricted vessel outline
try
    fig_constrictedbv = make_fig('roi_constricted_bv');
    fig_constrictedbv.update_figsize([8 6]);
    fig_constrictedbv.reset_axis();
    % Channel 1 background, green constricted, red dilated
    fig_constrictedbv.showrois(roilist, 1, ["manual_constricted_bv", "manual_dilated_bv"], ["-g", "-r"]);
    fig_constrictedbv.save2svg(save_dir);
catch ME
    warning(ME.identifier, 'Failed to generate roi_constricted_bv: %s', ME.message);
end

%% 5.1.3 PVS_dilated image + dilated and constricted vessel outline
try
    fig_dilatedpvs = make_fig('roi_dilated_pvs');
    fig_dilatedpvs.update_figsize([8 6]);
    fig_dilatedpvs.reset_axis();
    % Channel 2 background
    fig_dilatedpvs.showrois(roilist, 2, ["manual_dilated_pvs", "manual_constricted_pvs"], ["-r", "-g"]);
    fig_dilatedpvs.save2svg(save_dir);
catch ME
    warning(ME.identifier, 'Failed to generate roi_dilated_pvs: %s', ME.message);
end

%% 5.1.4 PVS_constricted image + dilated and constricted vessel outline
try
    fig_constrictedpvs = make_fig('roi_constricted_pvs');
    fig_constrictedpvs.update_figsize([8 6]);
    fig_constrictedpvs.reset_axis();
    % Channel 2 background
    fig_constrictedpvs.showrois(roilist, 2, ["manual_constricted_pvs", "manual_dilated_pvs"], ["-g", "-r"]);
    fig_constrictedpvs.save2svg(save_dir);
catch ME
    warning(ME.identifier, 'Failed to generate roi_constricted_pvs: %s', ME.message);
end

%% 5.2 Dip in, Dip out
try
    fig_dipinout = make_fig('roi_dipinout');
    fig_dipinout.update_figsize([8 6]);
    fig_dipinout.reset_axis();
    % Channel 1 background, green dipin, red dipout
    fig_dipinout.showrois(roilist, 1, ["manual_dipin", "manual_dipout"], ["-g", "-r"]);
    fig_dipinout.save2svg(save_dir);
catch ME
    warning(ME.identifier, 'Failed to generate roi_dipinout: %s', ME.message);
end

%% 5.3 PAX, Dilated, Constricted
try
    fig_pax = make_fig('roi_pax');
    fig_pax.update_figsize([8 6]);
    fig_pax.reset_axis();
    % Channel 3 (RGB)
    % Note: 'pax' is created by analysis_pax_fwhm (likely without 'manual_')
    % But dipin/out and bv definitions are manual
    fig_pax.showrois(roilist, 3, ["pax", "manual_dipin", "manual_dipout", "manual_constricted_bv", "manual_dilated_bv"], ["-c", "-y", "-y", "y", "y"]);
    fig_pax.save2svg(save_dir);
catch ME
    warning(ME.identifier, 'Failed to generate roi_pax: %s', ME.message);
end

end
