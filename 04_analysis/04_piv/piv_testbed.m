addpath(genpath(pwd));
sessiondir ='G:\tmp\01_igkltdt\hql086\250909_hql086_sleep\HQL086_sleep250909_005_piv';
session = ECSSession(sessiondir);
session = session.load_primary_results();
sleep_integrate = state_integration(sessiondir);
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
img_state = state_image(sleep_integrate);
img_state.get_state_indices(twophoton_processed.t_axis,twophoton_processed.outfps);
img_state.param.pixel2um = twophoton_processed.pixel2um;

%%
na_trans = get_stateframes(img_state.state_idx.na_trans,twophoton_processed.ch1);
ra_trans = get_stateframes(img_state.state_idx.ra_trans,twophoton_processed.ch1);
awake = get_stateframes(img_state.state_idx.awake,twophoton_processed.ch1);





%%
of = struct();
of.state = ra_trans{1};
of.pre = opticalflow_preprocess(of.state);
%% Check what will be feed in
of.input =of.pre(:,:,81:112);
sliceViewer(of.input);
%%
session.roilist.addormodifyroi(of.pre ,'piv','rectangle')
piv_mask = session.roilist.getvertices('piv');
of.uv = opticalflow_wlet(of.pre(:,:,60:150),"pyramid_levels",1,"smoothness",50);

%% Example of using the new ensemble multipass optical flow
[of.uv_ens, of.corr_ens] = corr_ensemble(of.pre(:,:,81:112), ...
    'window_sizes', [40 20], ...
    'repeat', 1, ...
    'do_pad', 1, ...
    'use_gpu', true);
%%
of.uv = opticalflow_wlet(of.pre(:,:,:), "roirect", roi_vertices2rect(piv_mask));
%%

pivensemble(of.pre(:,:,81:111))

%% select range while observing quivermap
of.range = opticalflow_viewer(of.state,of.uv,"block_size",11,"scale",10);
of.range = [81,111];

%%
figure()
imagesc(piv_divergence(of.uv_ens(:,:,:)*size(of.pre,3)*1000*twophoton_processed.pixel2um))
axis image
%%
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
figure()
imshow(of.pre(:,:,of.range(2)))
hold on
opticalflow_quiver(of.uv_ens,"block_size",5,"time_range",of.range,"scale",50,"linewidth",1)
%%
figure()
imshow(of.pre(:,:,of.range(2)))
hold on
piv_plotquiver(of.uv_ens,"block_size",1,"scale",size(of.pre,3)*2,"linewidth",1)