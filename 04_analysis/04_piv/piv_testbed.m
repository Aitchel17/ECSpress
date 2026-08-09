addpath(genpath(pwd));

%% ===== MANUAL TIFF MODE (no session): pick a TIFF, run PIV, view =====
% Quick PIV test on any TIFF stack, independent of ECSSession.
[tf, tp] = uigetfile({'*.tif;*.tiff', 'TIFF stack'}, 'Select a TIFF stack for PIV test');
tiffpath = fullfile(tp, tf);
% read multipage TIFF -> H x W x N double
tinfo  = imfinfo(tiffpath);
mstack = zeros(tinfo(1).Height, tinfo(1).Width, numel(tinfo));
tobj   = Tiff(tiffpath, 'r');
for k = 1:numel(tinfo); tobj.setDirectory(k); mstack(:,:,k) = double(tobj.read()); end
tobj.close();
% spatial (motion-preserving) preprocess + optional frame window
mpre = piv_preprocess(mstack);
mwin = 1:size(mpre,3);                       % e.g. 81:112 for a subrange
% ensemble PIV
m_cp = piv_corr_ensemble(mpre(:,:,mwin), ...
    window_sizes=[64 32; 32 16; 16 8], repeat=1, do_pad=1, use_gpu=false);
[m_cp.utable, m_cp.vtable] = piv_validate(m_cp.utable, m_cp.vtable, m_cp.corr);
[m_uv, m_corr] = vfield_stamp(m_cp);
%% divergence + quiver view
m_div = vfield_divergence_fd(m_uv);
figure; imagesc(m_div); axis image; colormap turbo; colorbar; title('divergence')
figure; imshow(mpre(:,:,mwin(end)), []); hold on
vfield_plotquiver(m_uv, block_size=1, scale=20, linewidth=1); title('PIV vectors')

%% ===== SESSION PIPELINE (below) =====
sessiondir ='E:\02_egfptdt\hql099\260601_hql99_sleep\HQL99_sleep260601_005';
session = ECSSession(sessiondir);
session = session.load_primary_results();
%todo 2026 04 28
% 1. state transition
% 2. PIV
% 3. Make offset video
% 1. Load data & 3. Load processed data (Integrated via ECSSession)
session.stackch1 = session.loadstack('ch1');
session.stackch2 = session.loadstack('ch2');
% 2. Twophoton data FPS matching & preprocessing
% Note: twophoton_preprocess expects a struct with stackch1/2 and img_param.
% ECSSession object works here as it has these properties.
twophoton_processed = twophoton_preprocess(session);
%%
sleep_integrate = state_integration(sessiondir);
img_state = state_image(sleep_integrate);
img_state.get_state_indices(twophoton_processed.t_axis,twophoton_processed.outfps);
img_state.param.pixel2um = twophoton_processed.pixel2um;

%%
na_trans = get_stateframes(img_state.state_idx.na_trans,twophoton_processed.ch1);
ra_trans = get_stateframes(img_state.state_idx.ra_trans,twophoton_processed.ch1);
awake = get_stateframes(img_state.state_idx.awake,twophoton_processed.ch1);
%%

nr_trans = get_stateframes(img_state.state_idx.nr_trans,twophoton_processed.ch1,150);



%%
of = struct();
%%
of.state = nr_trans{2};
of.pre = piv_preprocess(of.state);
figure()
sliceViewer(of.pre);
%% Check what will be feed in
of.range = [247 326]; % constriction [326,363], dilation [231 330] , nrem stable [20 120]
of.input =of.pre(:,:,of.range(1):of.range(2));

%%
figure()
sliceViewer(of.pre);


%% PIV masks (drawn on the preprocessed stack) -- applied right before ensemble
% piv_include : analyze ONLY inside this region  (true = keep)
% piv_exclude : drop these regions from analysis  (true = remove)
session.roilist.addormodifyroi(of.pre, 'piv_include', 'polygon');
session.roilist.addormodifyroi(of.pre, 'piv_exclude', 'polygon');
%%
incl = session.roilist.getmask('piv_include');
excl = session.roilist.getmask('piv_exclude');
if ~any(incl(:)); incl = true(size(incl)); end   % empty include -> whole frame
of.exclmask = ~incl | excl;                        % PIVlab convention: true = excluded

%% Example of using the new ensemble multipass optical flow
of.cp = piv_corr_ensemble(of.input, ...
    'window_sizes', [40 20], ...
    'repeat', 1, ...
    'do_pad', 1, ...
    'exclmask', of.exclmask, ...
    'use_gpu', true);
[of.cp.utable, of.cp.vtable] = piv_validate(of.cp.utable, of.cp.vtable, of.cp.corr);
[of.uv_ens, of.corr_ens] = vfield_stamp(of.cp);

%% ensemble piv plot
figure()
imshow(of.pre(:,:,of.range(1)))
hold on
vfield_plotquiver(of.uv_ens,"block_size",1,"time_range",of.range,"scale",10,"linewidth",1,...
    "head_len",5,"head_width",5,"filled_head",true ,"head_abs",true)
%% ensemble divergence
nfr   = size(of.input, 3);                                                       % window frame count
[of.tiles3, of.ti3] = sector_block(excl, 3*vfield_grid_step(of.uv_ens));
of.tiles3(of.exclmask) = 0;
of.tiles3 = sector_validate(of.tiles3, ~isnan(of.uv_ens(:,:,1)), of.ti3.n_labels, 9, 3);
m_div = sector_paint(of.tiles3, vfield_divergence(of.uv_ens, of.tiles3, of.ti3.n_labels)) * nfr;    % dimensionless (total areal strain)
figure; imagesc(m_div); axis image; colormap turbo; colorbar; title('divergence (block 2x2, dimensionless: \DeltaArea/Area)')

%% perivascular polar divergence  (rings from vessel wall x 4 quarters; min 4 vectors/cell)
nfr   = size(of.input, 3);                                        % window frame count
vmask   = session.roilist.getmask('piv_exclude');                 % vessel
uv_tot  = of.uv_ens * nfr;                                         % TOTAL displacement over window (px) -> all units consistent
mindisp = 0.00;                                                     % px: cell mean TOTAL |disp| below this -> NaN (~= prior 0.005 px/frame)
[pdiv, pcells, pinfo] = vfield_divergence_polar(uv_tot, vmask, ...
    'n_sectors', 8, 'min_valid', 4, 'exclmask', of.exclmask, 'min_disp', mindisp, ...
    'method', 'radial');   % v_r per (ring,sector) -> dv_r/dr across contours (no tangential)
%%
[pdiv, pcells, pinfo] = vfield_divergence_polar(uv_tot, vmask, ...
    'n_sectors', 8, 'min_valid', 4, 'min_disp', mindisp, ...
    'method', 'radial'); 

% pdiv is already dimensionless total areal strain (input field was total displacement)
bg = mat2gray(of.pre(:,:,of.range(1)));
figure; imshow(repmat(bg, [1 1 3])); hold on                      % truecolor bg (colormap-independent)
hd = imagesc(pdiv); set(hd, 'AlphaData', 0.55 * ~isnan(pdiv));
axis image; colormap turbo; colorbar
vfield_plotquiver(of.uv_ens, "block_size",1, "time_range",of.range, "scale",5, ...
    "filled_head",true, "head_abs",true, "head_len",5, "head_width",5, "linewidth",1)
%title(sprintf('perivascular divergence: %d rings x %d quarters (dimensionless)', pinfo.n_rings, pinfo.n_sectors))

%% perivascular radial velocity (v_r) map   -- in/out per cell (no tangential, no gradient)
bg = mat2gray(of.pre(:,:,of.range(1)));
figure; imshow(repmat(bg, [1 1 3])); hold on                       % truecolor bg
hv = imagesc(pinfo.vr_map); set(hv, 'AlphaData', 0.55 * ~isnan(pinfo.vr_map));
axis image; colormap turbo; colorbar
cmax = max(abs(pinfo.vr_map(:)), [], 'omitnan');                   % symmetric caxis, 0 = neutral
if ~isnan(cmax) && cmax > 0; clim([-cmax cmax]); end
vfield_plotquiver(of.uv_ens, "block_size",1, "time_range",of.range, "scale",5, ...
    "filled_head",true, "head_abs",true, "head_len",5, "head_width",5, "linewidth",1)
title('perivascular radial velocity v_r  (+ outward / - inward, total px)')

%% wavelet based method
session.roilist.addormodifyroi(of.pre ,'piv','rectangle')
piv_mask = session.roilist.getvertices('piv');
of.uv = opticalflow_wlet(of.pre(:,:,60:150),"pyramid_levels",1,"smoothness",50);

of.uv = opticalflow_wlet(of.pre(:,:,:), "roirect", roi_vertices2rect(piv_mask));
%%

pivensemble(of.pre(:,:,81:111))

%% select range while observing quivermap
of.range = vfield_viewer(of.state,of.uv,"block_size",11,"scale",10);
of.range = [51,90];

%%
figure()
imagesc(vfield_divergence_fd(of.uv_ens(:,:,:)*size(of.pre,3)*1000*twophoton_processed.pixel2um, 'exclmask', of.exclmask))
axis image
%%
figure()
imagesc(sum(squeeze(of.uv(:,:,of.range(1):of.range(2),2)),3))
axis image
%%
figure()
plot(of.uv(:,50,85,1))
%%
figure()
sliceViewer(cat(3,of.pre(:,:,1),of.pre(:,:,end)))
%%

%% divergence + quiver view  (2x2 block LSQ gradient; incomplete stencils -> NaN)
nfr   = size(of.input, 3);                                                       % window frame count
[of.tiles2, of.ti2] = sector_block(excl, 2*vfield_grid_step(of.uv_ens));
of.tiles2(of.exclmask) = 0;
of.tiles2 = sector_validate(of.tiles2, ~isnan(of.uv_ens(:,:,1)), of.ti2.n_labels, 4, 2);
m_div = sector_paint(of.tiles2, vfield_divergence(of.uv_ens, of.tiles2, of.ti2.n_labels)) * nfr;    % dimensionless (total areal strain)


figure; imagesc(m_div); axis image; colormap turbo; colorbar; title('divergence (block 2x2, dimensionless: \DeltaArea/Area)')
figure; imshow(mpre(:,:,mwin(end)), []); hold on
vfield_plotquiver(of.uv_ens, block_size=1, scale=20, linewidth=1); title('PIV vectors')

%%
figure()
imshow(of.pre(:,:,of.range(2)))
hold on
vfield_plotquiver(of.uv_ens,"block_size",1,"scale",size(of.pre,3)*2,"linewidth",1)