
bv_boundary = session.pax_fwhm.reconstruction(session.pax_fwhm.mask.upline+session.pax_fwhm.mask.downline);
pvs_boundary = session.pax_fwhm.reconstruction(session.pax_fwhm.mask.pvs_upline+session.pax_fwhm.mask.pvs_downline);
%%
channel_n = 4;
pax_fwhm_result = repmat(zeros(size(twophoton_processed.ch1)),[1,1,1,channel_n]);
pax_fwhm_result(:,:,:,1) = twophoton_processed.ch2;
pax_fwhm_result(:,:,:,2) = twophoton_processed.ch1;

%% blank generation
total_boundary = bv_boundary+pvs_boundary;
blank_idx = repmat(total_boundary>0, [1,1,1,channel_n]);
pax_fwhm_result(blank_idx) = 0;
%% Implant boundary
pax_fwhm_result(:,:,:,3) = bv_boundary*2048;
pax_fwhm_result(:,:,:,4) = pvs_boundary*2048;
%%
[d,b] = io_readtifftag()

%%
pax_fwhm_result = permute(pax_fwhm_result,[1,2,4,3]);
util_checkstack(pax_fwhm_result)
%%
io_postsavetiff(pax_fwhm_result,fullfile(sessiondir,"fwhm_tracking.tif"),[session.img_param.pixel2um,session.img_param.pixel2um,twophoton_processed.outfps])