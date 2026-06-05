% PIV good list
% 'G:\tmp\01_igkltdt\hql086\250909_hql086_sleep\HQL086_sleep250909_005'
% G:\tmp\01_igkltdt\hql086\250909_hql086_sleep\HQL086_sleep250909_007
addpath(genpath(pwd));
clc, clear
% for i = [2,3,4,5,6,7,8]
% sessiondir = strcat('G:\tmp\00_igkl\hql090\251020_hql090_whiskerb\HQL090_whiskerb251020_00',string(i));
% %%
sessiondir = 'G:\tmp\00_igkl\hql097\260321_hql097_sleep\HQL097_sleep260321_010';
%%
session = ECSSession(sessiondir);
session = session.load_primary_results();
state_integrate = state_integration(sessiondir);


%
paxfwhm_state = state_linefwhm(state_integrate);
paxfwhm_state.get_state_indices(session.pax_fwhm.t_axis,session.pax_fwhm.param.fs);

paxfwhm_state = state_linefwhm(state_integrate);
paxfwhm_state.get_state_indices(session.pax_fwhm.t_axis,session.pax_fwhm.param.fs);
contents_types = {["thickness","eps","bv","totalpvs","dynamic_pvs","static_pvs"],...
                ["displacement","dynamicpvs","staticpvs","dynamicbv","staticbv"]};
%%


for typeidx = 1:numel(contents_types)
    contents_type = contents_types{typeidx};
    d_type = contents_type(1);
    for contents_idx =  2:numel(contents_type)
        content_name = contents_type(contents_idx);
        name = strcat(d_type,'_',content_name);
        fprintf('Processing %s...\n', name);
        data = session.pax_fwhm.(d_type).(content_name);
        paxfwhm_state.get_summary(name, data);
        paxfwhm_state.get_powerdensity(name, data);
        paxfwhm_state.decompose_signal(name, data);
        paxfwhm_state.get_pppt_decomposition(name);
        paxfwhm_state.get_transitionsummary(name, data);
    end
end
% end
paxfwhm_state.save2disk('paxfwhm_state',state_integrate.dir_struct.stateanalysis)

%% Image
analog = load(state_integrate.dir_struct.analog);
analog = analog.primary_analog;
%todo 2026 04 28
% 1. state transition
% 2. PIV
% 3. Make offset video
% 1. Load data & 3. Load processed data (Integrated via ECSSession)
session = ECSSession(sessiondir);
session = session.load_primary_results();
session.stackch1 = session.loadstack('ch1');
session.stackch2 = session.loadstack('ch2');
% 2. Twophoton data FPS matching & preprocessing
% Note: twophoton_preprocess expects a struct with stackch1/2 and img_param.
% ECSSession object works here as it has these properties.
twophoton_processed = twophoton_preprocess(session);
%%
img_state = state_image(state_integrate);
img_state.get_state_indices(twophoton_processed.t_axis,twophoton_processed.outfps);
img_state.param.pixel2um = twophoton_processed.pixel2um;
%session.polarcluster = analysis_clusterpolar(session.pax_fwhm, twophoton_processed, session.dir_struct.primary_analysis);
%%
na_trans = get_stateframes(img_state.state_idx.na_trans,twophoton_processed.ch1);
ra_trans = get_stateframes(img_state.state_idx.ra_trans,twophoton_processed.ch1);

%%
result = opticalflow_preprocess(ra_trans{1});
%%
util_checkstack(result)
%%
of_data =opticalflow_wlet(result);
%%
opticalflow_quiver(of_data(:,:,1,:),"block_size",5,"scale",10)
%%

op_range = opticalflow_viewer(result,of_data,"block_size",5,scale=10);
%%
figure()
imshow(result(:,:,op_range(2)))

opticalflow_quiver(of_data,"block_size",5,"time_range",op_range,"scale",1,"linewidth",1)

%%
I0 = (result(:,:,op_range(1))-min(result(:,:,op_range(1))))./max(result(:,:,op_range(1)));
I1 = (result(:,:,op_range(2))-min(result(:,:,op_range(2))))./max(result(:,:,op_range(2)));
%%

D = I1 - I0;

figure;
imagesc(log(D-min(D,[],"all")));
axis image
%%
logd = log(D-min(D,[],"all"));
%%
figure()
imshow(cat(3,I0,logd*10,I0))


%%
figure()
imagesc(norm_result1-norm_result2)
%%

rgb_of = cat(3,result(:,:,op_range(2)),result(:,:,op_range(1)),zeros(size(result(:,:,op_range(2)))));

%%
sliceViewer(cat(3,result(:,:,op_range(1)),result(:,:,op_range(2))))
%%
session.roilist.add 


optical
%%
figure()
edge(result(:,:,op_range(1)),"canny")
%%
figure()
edge(result(:,:,op_range(2)),"canny")
%%
figure()
imshow(cat(3,edge(result(:,:,op_range(2)),"canny"),edge(result(:,:,op_range(1)),"canny"),zeros(size(result(:,:,op_range(2))))))


%%
trans_img.dir = fullfile(session.dir_struct.primary_analysis,"trans");
mkdir(trans_img.dir)
io_postsavetiff(ra_trans{1},fullfile(session.dir_struct.primary_analysis,'ratrans.tif'),...
    [twophoton_processed.pixel2um,twophoton_processed.pixel2um,twophoton_processed.fps])
%%
bv = paxfwhm_state.transition.thickness_bv.data{:,1};
%%
figure()
plot(bv)

%%
bv = paxfwhm_state.transition.thickness_bv.data{1};
%%
cla
%%
figure()

plot(session.pax_fwhm.t_axis,session.pax_fwhm.thickness.bv)




%% debugging script


figure()
%%
cla
plot(analog.ball.rs_taxis,analog.ball.resampled_ball)
%%

plot(session.pax_fwhm.t_axis,session.pax_fwhm.thickness.bv)

hold on
%%
plot(session.pax_fwhm.t_axis,session.pax_fwhm.thickness.pvschanges_total,Color='g')


%%
move = state_integrate.time_table.;
y_pos = max(session.pax_fwhm.thickness.bv);
            for idx = 1:numel(move)
                stim_dur = move(idx);

                % Retrieve the Nx2 table/matrix of [StartTime, EndTime]
                time_intervals = move;

                % Fast plotting strategy: Create vectors separated by NaNs 
                % so we only call plot() once per duration type.
                num_intervals = size(time_intervals, 1);
                x_data = NaN(3, num_intervals);
                x_data(1, :) = time_intervals(:, 1)'; % StartTime
                x_data(2, :) = time_intervals(:, 2)'; % EndTime

                y_data = NaN(3, num_intervals);
                y_data(1:2, :) = y_pos; % Set the line height to 2.5

                % Plot all horizontal bars for this duration type
                plot(x_data(:), y_data(:), 'Color', 'g', ...
                     'LineWidth', 4);
            end

