sessiondir = 'G:\tmp\00_igkl\hql090\251016_hql090_sleep\HQL090_sleep251016_008';
session = ECSSession(sessiondir);
session = session.load_primary_results();
sleep_integrate = state_integration(sessiondir);

%% Paxfwhm
paxfwhm_state = state_linefwhm(sleep_integrate);
paxfwhm_state.get_state_indices(session.pax_fwhm.t_axis,session.pax_fwhm.param.fs);

contents_types = {["thickness","eps","bv","totalpvs","dynamic_pvs","static_pvs"],...
                ["displacement","dynamicpvs","staticpvs","dynamicbv","staticbv"]};

for typeidx = 1:numel(contents_types)
    contents_type = contents_types{typeidx};
    d_type = contents_type(1);
    fprintf('Processing %s...\n', name);
    for contents_idx =  2:numel(contents_type)
        content_name = contents_type(contents_idx);
        name = strcat(d_type,'_',content_name);
        data = session.pax_fwhm.(d_type).(content_name);
        paxfwhm_state.get_summary(name, data);
        paxfwhm_state.get_powerdensity(name, data);
        paxfwhm_state.decompose_signal(name, data);
        paxfwhm_state.get_pppt_decomposition(name);
        paxfwhm_state.get_transitionsummary(name, data);
    end
end
paxfwhm_state.save2disk('paxfwhm_sleep',sleep_integrate.dir_struct.stateanalysis)






%% Plot transition analysis examples
clee = color_lee;
%%
figure()
hold on

plot_state = "na_trans";
plot_target = ["thickness_bv", "thickness_totalpvs","thickness_eps"];
color_map = {clee.clist.red, clee.clist.green, clee.clist.orange};

for target_idx = 1:numel(plot_target)
    transition_table = paxfwhm_state.get_filtered_table('transition',plot_target(target_idx), plot_state);
    plot(mean(cell2mat(transition_table.data),1)/median(cell2mat(transition_table.data),"all"),"Color",color_map{target_idx});
end

%%

figure()
hold on

plot_state = "ra_trans";
plot_target = ["thickness_bv", "thickness_totalpvs"];
color_map = {clee.clist.red, clee.clist.coloredgreen};

for target_idx = 1:numel(plot_target)
    transition_table = paxfwhm_state.get_filtered_table('transition',plot_target(target_idx), plot_state);
    mean(cell2mat(transition_table.data),1);
    plot(mean(cell2mat(transition_table.data),1)/median(cell2mat(transition_table.data),"all"),"Color",color_map{target_idx});
end




%%
figure()
    transition_table1 = paxfwhm_state.get_filtered_table('transition',"thickness_ecschanges_residual", "nr_trans");
    data1 = mean(cell2mat(transition_table1.data),1);
    plot(data1, "Color",'r');





%%
figure()
    transition_table1 = paxfwhm_state.get_filtered_table('transition',"thickness_bv", "nr_trans");
    data1 = mean(cell2mat(transition_table1.data),1);

    transition_table2 = paxfwhm_state.get_filtered_table('transition',"thickness_totalpvs", "nr_trans");
    data2 = mean(cell2mat(transition_table2.data),1);
    plot(data1-data2,"Color",'r');



%%







%%
paxfwhm_state.save2disk(paxfwhm_state,sleep_integrate)


%%
bv_nrem_table = paxfwhm_state.get_filtered_table('state_summary','BV_thickness',"nrem");
bv_awake_table = paxfwhm_state.get_filtered_table('state_summary','BV_thickness',"awake");
bv_rem_table = paxfwhm_state.get_filtered_table('state_summary','BV_thickness',"rem");

% Define column to plot
target_col = 'raw_median'; % Variable to compare
% Extract data vectors
data_nrem = bv_nrem_table.(target_col);
data_awake = bv_awake_table.(target_col);
data_rem = bv_rem_table.(target_col);

% Plot
plot_bar_points({ data_awake, data_nrem, data_rem}, { 'Awake', 'NREM', 'REM'}, ...
    'YLabel', target_col, 'Title', ['Comparison: ' target_col], ...
    'Colors', [0 0 1; 0 1 0; 1 0 0]); % Blue, Red
%% Define column to plot
target_col = 'raw_median'; % Variable to compare
% Extract data vectors
pvs_nrem_table = paxfwhm_state.get_filtered_table('state_summary','PVStotal_thickness',"nrem");
pvs_awake_table = paxfwhm_state.get_filtered_table('state_summary','PVStotal_thickness',"awake");
pvs_rem_table = paxfwhm_state.get_filtered_table('state_summary','PVStotal_thickness',"rem");



data_pvsnrem = pvs_nrem_table.(target_col);
data_pvsawake = pvs_awake_table.(target_col);
data_pvsrem = pvs_rem_table.(target_col);

% Plot
plot_bar_points({ data_pvsawake, data_pvsnrem, data_pvsrem}, { 'Awake', 'NREM', 'REM'}, ...
    'YLabel', target_col, 'Title', ['Comparison: ' target_col], ...
    'Colors', [0 0 1; 0 1 0; 1 0 0]); % Blue, Red

%%
paxfwhm_state.band_decomposition






%%
% todo: 20260205
% 0. Change the prpoerty name... maybe summaryanalysis -> state_summary,
% poweranalysis powerdensitym, decomposition--> band_decomposition analysis
% transition
%  perhaps what happens is vessel dilation and ECS change is lagging
%  therefore summary would not reflect what actually happened during this
%  decomposition analysis is very necessary, the simple stat of arousal
%  state is meaningless as sleep is very fragmented and often previous
%  state is reflected to next state even for the long fluctuation

% 2. decomposition
%   butter filter artifact issue medianmovwindow?
%   might be better to make decompose_signal func to get bandinput --> no
%   its better not to use function to hide band parameter..

% 1. Add start time, end time, and composition to the summary analysis (Done)
% 2. Embed decomposition analysis(pp pt) to stateanalysis_1d (In progress)
% 2.1 Decompse signal (Done)
% 2.2 PP PT analysis (Done)

% 3. Transition analysis need to be seperated from summary_analysis and (Done)
% treated seperately



% 4. Saving mechanism (Done)
% 5. Integration (manually in script to show something tomorrow)
% 6. Embed plot_sleep_patches(gca, sleep_score); to make_fig

%% Todo: 20260207
% 0. plot script
% 1. integration script


%%

paxfwhm_state.get_powerdensity('bv_thickness',session.pax_fwhm.thickness.bv) % 0.04 Hz for 25 second window

%%



tmp.decomposed_bv = decompose_signal(session.pax_fwhm.thickness.bv,session.pax_fwhm.param.fs);
tmp.decomposed_pvs = decompose_signal(session.pax_fwhm.thickness.totalpvs,session.pax_fwhm.param.fs);

%%
figure()

plot(session.pax_fwhm.t_axis,tmp.decomposed_bv.continuous,'k');
hold on
plot(session.pax_fwhm.t_axis,tmp.decomposed_pvs.continuous ,'b');
%%
figure()

plot(session.pax_fwhm.t_axis,tmp.decomposed_bv.isf,'k');
hold on
plot(session.pax_fwhm.t_axis,tmp.decomposed_pvs.isf ,'b');



%%



tmp.data = tmp.decomposed_bv.vlf;
plot(session.pax_fwhm.t_axis,tmp.data,'k');


%%
figure()
plot(session.pax_fwhm.idx.lowerBVboundary - session.pax_fwhm.idx.upperBVboundary)
hold on
plot(session.pax_fwhm.thickness.bv)

%%
plot_sleep_patches(gca, sleep_score);






%%
data1d.thickness = session.pax_fwhm.thickness;
data1d.displacement = session.pax_fwhm.displacement;
data1d.idx = session.pax_fwhm.idx;
%%


StateAnalysis_1D(sleep_integrate,session)
%%
clee = color_lee;
%%
session.stackch1 = session.loadstack('ch1');
session.stackch2 = session.loadstack('ch2');
%% Two photon analysis time axis calculation
session.pax_fwhm.t_axis = linspace(session.img_param.imgstarttime,session.img_param.imgendtime,numel(session.pax_fwhm.idx.upperBVboundary));
sleep_score = load(fullfile(sessiondir,"peripheral","sleep_score.mat"));
%%



%%
% Calculate
% Find indices for all REM epochs
% Calculate State Analysis (Indices, Bouts, Sliced Data)
rem.pax_fwhm = add_state_analysis(session.pax_fwhm, sleep_score.REMTimes);
nrem.pax_fwhm = add_state_analysis(session.pax_fwhm, sleep_score.NREMTimes, sleep_score.uArousalTimes);
%% Big Chunk Spectral Analysis (300s)
param_fs = session.pax_fwhm.param.fs; % Moved definition of param_fs here to ensure it's available
rem.bigchunk = get_bigchunk_analysis(session.pax_fwhm, sleep_score, "REM", param_fs, 300);
nrem.bigchunk = get_bigchunk_analysis(session.pax_fwhm, sleep_score, "NREM", param_fs, 300);
awake.bigchunk = get_bigchunk_analysis(session.pax_fwhm, sleep_score, "Awake", param_fs, 300);
drowsy.bigchunk = get_bigchunk_analysis(session.pax_fwhm, sleep_score, "Drowsy", param_fs, 300);
awake.pax_fwhm = add_state_analysis(session.pax_fwhm, sleep_score.AwakeTimes);
drowsy.pax_fwhm = add_state_analysis(session.pax_fwhm, sleep_score.DrowsyTimes);
rem.img.param = session.img_param;
rem.img.stateidx = statebin2frame(sleep_score.REMTimes, session.img_param.taxis);
nrem.img.param = session.img_param;
nrem.img.stateidx = statebin2frame(sleep_score.NREMTimes, sleep_score.uArousalTimes, session.img_param.taxis);
awake.img.param = session.img_param;
awake.img.stateidx = statebin2frame(sleep_score.AwakeTimes, session.img_param.taxis);
drowsy.img.param = session.img_param;
drowsy.img.stateidx = statebin2frame(sleep_score.DrowsyTimes, session.img_param.taxis);
% bouts loc of image
nrem.img.bouts = stateidx2bouts(nrem.img.stateidx);
awake.img.bouts = stateidx2bouts(awake.img.stateidx);
rem.img.bouts = stateidx2bouts(rem.img.stateidx);
drowsy.img.bouts = stateidx2bouts(drowsy.img.stateidx);
% Spectrogram Calculation (Batch Process All 1D Data)
param_fs = session.pax_fwhm.param.fs;
rem.spectral_analysis = get_all_spectra(rem.pax_fwhm.bouts, session.pax_fwhm, param_fs);
nrem.spectral_analysis = get_all_spectra(nrem.pax_fwhm.bouts, session.pax_fwhm, param_fs);
awake.spectral_analysis = get_all_spectra(awake.pax_fwhm.bouts, session.pax_fwhm, param_fs);
drowsy.spectral_analysis = get_all_spectra(drowsy.pax_fwhm.bouts, session.pax_fwhm, param_fs);

%%
transition.roughawake = statebin2timetable(sleep_score.AwakeTimes, sleep_score.DrowsyTimes);
transition.roughnrem = statebin2timetable(sleep_score.NREMTimes, sleep_score.uArousalTimes);
transition.rem = statebin2timetable(sleep_score.REMTimes);
%%

%%
transition.window = 25;
transition.NA.timetable = get_transition(transition.window,transition.roughnrem, transition.roughawake);
transition.NA.fwhmidx = timearr2frame(session.pax_fwhm.t_axis, transition.NA.timetable);
transition.AN.timetable = get_transition(transition.window,transition.roughawake, transition.roughnrem);
transition.AN.fwhmidx = timearr2frame(session.pax_fwhm.t_axis, transition.AN.timetable);
transition.NR.timetable = get_transition(transition.window,transition.roughnrem, transition.rem);
transition.NR.fwhmidx = timearr2frame(session.pax_fwhm.t_axis, transition.NR.timetable);
transition.RA.timetable = get_transition(transition.window,transition.rem, transition.roughawake);
transition.RA.fwhmidx = timearr2frame(session.pax_fwhm.t_axis, transition.RA.timetable);
%%
tmp.dataname = 'ecschanges_residual'; % totalpvs, bv
tmp.transitname = 'RA';

figure('Name',strcat(tmp.dataname,'_',tmp.transitname))


tmp.sumarray = [];
cla
hold on
tmp.data = sgolayfilt(session.pax_fwhm.thickness.(tmp.dataname),3,5);
tmp.ndata = size(transition.(tmp.transitname).fwhmidx,1);
tmp.idx = transition.(tmp.transitname).fwhmidx;

if strcmp(tmp.dataname,'bv')
    sub_color = [0.9 0.8 0.8];
    main_color = 'r';
else
    sub_color = [0.8 0.9 0.8];
    main_color = 'g';
end

for idx = 1:tmp.ndata
    tmp.start = tmp.idx(idx,1);
    tmp.end = tmp.idx(idx,2);
    tmp.plottarget = tmp.data(1,tmp.start:tmp.end);
    tmp.sumarray = [tmp.sumarray; tmp.plottarget];
    tmp.plottaxis = linspace(-transition.window,transition.window,length(tmp.plottarget));
    plot(tmp.plottaxis,tmp.plottarget, color=sub_color)
end

tmp.meantarget = mean(tmp.sumarray,1);
tmp.sem = std(tmp.sumarray, 0, 1, 'omitnan') ./ sqrt(tmp.ndata);
tmp.t_score = tinv(0.975, tmp.ndata-1); % 95% CI
tmp.ci = tmp.t_score * tmp.sem;

tmp.upci = tmp.meantarget + tmp.ci;
tmp.downci = tmp.meantarget - tmp.ci;

patch([tmp.plottaxis, fliplr(tmp.plottaxis)], [tmp.upci, fliplr(tmp.downci)],main_color, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
plot(tmp.plottaxis,tmp.meantarget, main_color)
xlim([-transition.window transition.window])
xline(0)
%%



%%
transition.RA
transition.AN
transition.NR
%%





%%
if ~isfolder(fullfile(sessiondir,'sleep_analysis'))
    mkdir(fullfile(sessiondir,'sleep_analysis'))
end

%% Image frame location (kept separate for now as it uses img_param.taxis)

%% OVERVIEW


figure()
%%
cla
plot(session.pax_fwhm.t_axis,sgolayfilt(session.pax_fwhm.thickness.bvchanges,2,5),'Color',clee.clist.red)
hold on
plot(session.pax_fwhm.t_axis,sgolayfilt(session.pax_fwhm.thickness.pvschanges_total,2,5),'Color',clee.clist.darkgreen)

yl = ylim;
% Draw Patches for Sleep States
% Draw Patches for Sleep States
plot_sleep_patches(gca, sleep_score);







%%
focus = awake;
bouts = 5;
figure(name='bvchanges_nrem')
plot(focus.pax_fwhm.thickness(bouts).bvchanges,'color',[0.5 0.5 0.5])
hold on
plot(focus.spectral_analysis.thickness.bvchanges.decomposed.continuous.data{bouts},'k')
plot(focus.spectral_analysis.thickness.bvchanges.decomposed.vlf.data{bouts},'r')
plot(focus.spectral_analysis.thickness.bvchanges.decomposed.lf.data{bouts},'g')
%%
xlim([0 150])
ylim([-5 5])
%%
%%
figure(name='pvschanges_nrem')
plot(focus.pax_fwhm.thickness(bouts).pvschanges_total,'color',[0.5 0.5 0.5])
hold on
plot(focus.spectral_analysis.thickness.pvschanges_total.decomposed.continuous.data{bouts},'k')

plot(focus.spectral_analysis.thickness.pvschanges_total.decomposed.vlf.data{bouts},'r')
plot(focus.spectral_analysis.thickness.pvschanges_total.decomposed.lf.data{bouts},'g')
%%

figure()
plot(focus.pax_fwhm.thickness(bouts).ecschanges_residual,'color',[0.5 0.5 0.5])
hold on
plot(focus.spectral_analysis.thickness.ecschanges_residual.decomposed.continuous{bouts},'k')

plot(focus.spectral_analysis.thickness.ecschanges_residual.decomposed.vlf{bouts},'r')
plot(focus.spectral_analysis.thickness.ecschanges_residual.decomposed.lf{bouts},'g')


%%
figure()
plot(focus.pax_fwhm.displacement(1).dynamicpvs,'color',[0.5 0.5 0.5])
hold on
plot(focus.spectral_analysis.displacement.dynamicpvs.decomposed.continuous{1},'k')

plot(focus.spectral_analysis.displacement.dynamicpvs.decomposed.vlf{1},'r')
plot(focus.spectral_analysis.displacement.dynamicpvs.decomposed.lf{1},'g')

%%
xlim([0 150])
ylim([-5 5])
%%

figure(name='vlf bv pvs ecs awake')
plot(focus.spectral_analysis.thickness.bvchanges.decomposed.vlf.data{1},'r')

hold on

plot(focus.spectral_analysis.thickness.pvschanges_total.decomposed.vlf.data{1},'g')

plot(focus.spectral_analysis.thickness.ecschanges_residual.decomposed.vlf.data{1},'color',[0.5 0.5 0.5])


ylim([-4 4])
xlim([0 150])
%%



mean(nrem.spectral_analysis.thickness.bv.peak_analysis.lf.p2t_avg)
mean(awake.spectral_analysis.thickness.bv.peak_analysis.lf.p2t_avg)


%%
clc
mean(awake.spectral_analysis.thickness.bv.peak_analysis.continuous.p2t_avg,'omitmissing')
mean(drowsy.spectral_analysis.thickness.bv.peak_analysis.continuous.p2t_avg,'omitmissing')
mean(nrem.spectral_analysis.thickness.bv.peak_analysis.continuous.p2t_avg,'omitmissing')
mean(rem.spectral_analysis.thickness.bv.peak_analysis.continuous.p2t_avg,'omitmissing')
%%
clc
mean(awake.spectral_analysis.thickness.bv.peak_analysis.vlf.p2t_avg,'omitmissing')
mean(drowsy.spectral_analysis.thickness.bv.peak_analysis.vlf.p2t_avg,'omitmissing')
mean(nrem.spectral_analysis.thickness.bv.peak_analysis.vlf.p2t_avg,'omitmissing')
mean(rem.spectral_analysis.thickness.bv.peak_analysis.vlf.p2t_avg,'omitmissing')
%%
mean(awake.spectral_analysis.thickness.bv.peak_analysis.lf.p2t_avg,'omitmissing')
mean(drowsy.spectral_analysis.thickness.bv.peak_analysis.lf.p2t_avg,'omitmissing')
mean(nrem.spectral_analysis.thickness.bv.peak_analysis.lf.p2t_avg,'omitmissing')
mean(rem.spectral_analysis.thickness.bv.peak_analysis.lf.p2t_avg,'omitmissing')
%%

mean(awake.spectral_analysis.thickness.bv.peak_analysis.continuous.p2t_avg,'omitmissing')
mean(drowsy.spectral_analysis.thickness.bv.peak_analysis.continuous.p2t_avg,'omitmissing')
mean(nrem.spectral_analysis.thickness.bv.peak_analysis.continuous.p2t_avg,'omitmissing')
mean(rem.spectral_analysis.thickness.bv.peak_analysis.continuous.p2t_avg,'omitmissing')

%%
clc
mean(awake.spectral_analysis.thickness.pvschanges_total.peak_analysis.lf.p2t_avg,'omitmissing')
mean(drowsy.spectral_analysis.thickness.pvschanges_total.peak_analysis.lf.p2t_avg,'omitmissing')
mean(nrem.spectral_analysis.thickness.pvschanges_total.peak_analysis.lf.p2t_avg,'omitmissing')
mean(rem.spectral_analysis.thickness.pvschanges_total.peak_analysis.lf.p2t_avg,'omitmissing')
%%
clc
mean(awake.spectral_analysis.thickness.pvschanges_total.peak_analysis.vlf.p2t_avg,'omitmissing')
mean(drowsy.spectral_analysis.thickness.pvschanges_total.peak_analysis.vlf.p2t_avg,'omitmissing')
mean(nrem.spectral_analysis.thickness.pvschanges_total.peak_analysis.vlf.p2t_avg,'omitmissing')
mean(rem.spectral_analysis.thickness.pvschanges_total.peak_analysis.vlf.p2t_avg,'omitmissing')


%%

%%
focus = awake;

[c, lags] = xcorr(focus.spectral_analysis.thickness.bvchanges.decomposed.vlf{1}, focus.spectral_analysis.thickness.pvschanges_total.decomposed.vlf{1});
%%
figure()
%%
plot(lags, c)



%%
%% Plotting Summary
figure()
hold on

% Define target for plotting (Example: Awake BV Changes)
% Access path: state.spectral_analysis.structName.fieldName.summary
target_summary = awake.spectral_analysis.thickness.bvchanges.spec_summary;
target_list = awake.spectral_analysis.thickness.bvchanges.list;

% Plot individual (log-log)
for i = 1:numel(target_list)
    loglog(target_list(i).F, target_list(i).S, 'Color', [0.8 0.6 0.6]) % light red
end
% Plot average
loglog(target_summary.F, target_summary.S, 'r', 'LineWidth', 2)

xlabel('Frequency (Hz)'); ylabel('Power');
title('Average spectral density (Awake - BV)');


% Plot individual (log-log) - NREM
target_summary_nrem = nrem.spectral_analysis.thickness.bvchanges.spec_summary;
target_list_nrem = nrem.spectral_analysis.thickness.bvchanges.list;

for i = 1:numel(target_list_nrem)
    loglog(target_list_nrem(i).F, target_list_nrem(i).S, 'Color', [0.6 0.8 0.6]) % light green
end
% Plot average
loglog(target_summary_nrem.F, target_summary_nrem.S, 'g', 'LineWidth', 2)
set(gca, 'XScale', 'linear', 'YScale', 'linear')
% Note: user previously requested linear scale here, but loglog implies log.
% Resetting to log since mostly spectral plots are log-log.
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('Frequency (Hz)'); ylabel('Power');
title('Average spectral density (Awake vs NREM BV)');

%% mtspectrum based analysis of fwhm thickness

%% spectral analysis

focus = drowsy.fwhmloc_bouts;
cnt = 1;
for i = 1:numel(focus)
    loc = focus{i};
    if numel(loc)>eventlen_thr
        disp(i)
        % Preprocessing: Sgolay filter + Detrend
        bv_signal = session.pax_fwhm.thickness.bvchanges(loc);
        % bv_signal = detrend(sgolayfilt(bv_signal, 3, 5));
        bv_spectrogram = get_spectrogram(session.pax_fwhm.param.fs,bv_signal);
        bv_spectrogram.boutidx = i;
        bv_spectrogramlist(cnt) = bv_spectrogram;
        cnt = cnt + 1;
    end
end


%%
figure(name='nrem_pvs')
hold on
clc
for i = 1:numel(focus)
    loc = focus{i};
    if numel(loc)>eventlen_thr
        disp(i)
        % Preprocessing: Sgolay filter + Detrend
        pvs_signal = session.pax_fwhm.thickness.pvschanges_total(loc);
        % bv_signal = detrend(sgolayfilt(bv_signal, 3, 5));
        bv_spectrogram = get_spectrogram(session.pax_fwhm.param.fs,pvs_signal);
        loglog(bv_spectrogram.F,bv_spectrogram.S)
    end
end

%%

% Calculate PSD for each NREM bout
nrem_psd = struct('bv', {}, 'pvs', {});
fs = session.pax_fwhm.param.fs; % Calculate fs from actual time axis

for i = 1:numel(fwhmloc.nrem_bouts)
    loc = fwhmloc.nrem_bouts{i};

    % Preprocessing: Sgolay filter + Detrend
    bv_signal = session.pax_fwhm.thickness.bvchanges(loc);
    bv_signal = detrend(sgolayfilt(bv_signal, 3, 11));

    pvs_signal = session.pax_fwhm.thickness.pvschanges_dynamic(loc);
    pvs_signal = detrend(sgolayfilt(pvs_signal, 3, 11));

    % BV power spectrum
    [pxx_bv, f] = pwelch(bv_signal-mean(bv_signal), [], [], [], fs);
    nrem_psd(i).bv.Pxx = pxx_bv;
    nrem_psd(i).bv.F = f;

    % PVS power spectrum
    [pxx_pvs, f] = pwelch(pvs_signal-mean(pvs_signal), [], [], [], fs);
    nrem_psd(i).pvs.Pxx = pxx_pvs;
    nrem_psd(i).pvs.F = f;
end


%%

fig = figure('Name','Thickness plot');
t = tiledlayout(fig,3,1);
loc = fwhmloc.nrem;
ax.idx = nexttile(t, 1);
hold on
imagesc(ax.idx, t_axis, 1:size(session.pax_fwhm.kymograph.kgph_pvs_processed,1), session.pax_fwhm.kymograph.kgph_pvs_processed)
plot(ax.idx,t_axis,session.pax_fwhm.idx.clean_lowerBVboundary,'r')
plot(ax.idx,t_axis,session.pax_fwhm.idx.clean_pvsdownedge_idx,'b')
plot(ax.idx,t_axis,session.pax_fwhm.idx.clean_upperBVboundary,'r')
plot(ax.idx,t_axis,session.pax_fwhm.idx.clean_pvsupedge_idx,'b')


ax.upidx = nexttile(t, 3);
plot(ax.upidx,t_axis(loc),session.pax_fwhm.idx.clean_upperBVboundary(loc),'r')
hold on
plot(ax.upidx,t_axis(loc),session.pax_fwhm.idx.clean_pvsdownedge_idx(loc),'b')

ax.thickness = nexttile(t, 2);
plot(ax.thickness,t_axis(loc),session.pax_fwhm.thickness.bvchanges(loc),'r')
hold on
plot(ax.thickness,t_axis(loc),session.pax_fwhm.thickness.pvschanges_dynamic(loc),'b')

all_axes = struct2cell(ax);
all_axes = [all_axes{:}]; % Convert cell array to object array
linkaxes(all_axes, 'x');



%%

channel = ["bv","csf"];
state = ["rem","nrem","awake"];

%%
fig = figure();
fig.Units = 'inches';
fig.Position = [21 1 18 7]; % [left bottom width height]
t = tiledlayout(numel(channel),numel(state), 'TileSpacing', 'none', 'Padding', 'none');

img_handles = gobjects(numel(channel),numel(state));
%
ax = struct();
for ch =1:numel(channel)
    for st = 1:numel(state)
        fieldname =  strcat(channel(ch),"_",state(st));
        % Use nexttile to get the next tile in the layout
        ax.(fieldname) = nexttile(t);
        current_locs = imgloc.(state(st));
        if isempty(current_locs)
            title(ax.(fieldname), [fieldname " (Empty)"], 'Interpreter', 'none');
            continue;
        end

        if ch == 1
            img_data = mean(session.stackch1(:,:,current_locs), 3);
        elseif ch == 2
            img_data = mean(session.stackch2(:,:,current_locs), 3);
        end

        % Contrast adjustment: 10th to 90th percentile
        clim_val = prctile(img_data(:), [10 99]);
        img_handles(ch, st) = imagesc(ax.(fieldname), img_data, clim_val);
        axis(ax.(fieldname), 'image'); % Fix aspect ratio
        title(ax.(fieldname), fieldname, 'Interpreter', 'none');
    end
end
%% Sync all axes for simultaneous panning/zooming
all_axes = struct2cell(ax);
all_axes = [all_axes{:}]; % Convert cell array to object array
linkaxes(all_axes, 'xy');
%%
%% Interactive Marker Control (Click on Figure)
fprintf('--- Interactive Marker Mode ---\n');
fprintf('Click on any image to mark position on all axes.\n');
fprintf('Press Enter (key) to stop.\n');

marker_handles = [];
while true
    try
        [x_val, y_val, button] = ginput(1);
    catch
        break; % Handle figure closure
    end

    % Break if Enter (empty button) or any key press (button > 3 typically)
    % Mouse clicks are 1 (Left), 2 (Middle), 3 (Right)
    if isempty(button) || ~ismember(button, [1, 2, 3])
        fprintf('Exiting marker mode.\n');
        break;
    end

    % Remove old markers
    delete(marker_handles);
    marker_handles = gobjects(0);

    % Add new markers to all axes
    for i = 1:numel(all_axes)
        try
            curr_ax = all_axes(i);
            hold(curr_ax, 'on');
            h = plot(curr_ax, x_val, y_val, 'rx', 'MarkerSize', 8, 'LineWidth', 1);
            marker_handles(end+1) = h;
        catch
            % Handle closed axes
        end
    end
    drawnow;
end

%%

imagesc(ans)
axis image
%%




%%
session.pax_fwhm.kymograph.kgph_lumen_processed(:,drowsy_locs)


%%
clee = color_lee;
%%
thickness_fig = make_fig('Thickness figure');
%%
thickness_fig.hold_axis(true)
thickness_fig.resolution = session.img_param.pixel2um;
thickness_fig.update_figsize([15,3])
thickness_fig.monitor_xyinch = [21, 5];
thickness_fig.reset_axis
thickness_fig.loc.x = t_axis;
xlim(thickness_fig.ax,[0 max(t_axis)])
thickness_fig.plot_line(session.pax_fwhm.thickness.bvchanges,'r')
thickness_fig.plot_line(session.pax_fwhm.thickness.pvschanges_total, clee.clist.darkgreen)
thickness_fig.put_xaxistitle('Time (sec)')
thickness_fig.put_yaxistitle('Thickness change (um)')
%%
kymograph_fig = make_fig('Kymograph figure');
kymograph_fig.resolution = session.img_param.pixel2um;
kymograph_fig.update_figsize([15,3])
kymograph_fig.monitor_xyinch = [21, 5];
kymograph_fig.reset_axis
kymograph_fig.loc.x = t_axis;

kymograph_fig.plot_kymograph(session.pax_fwhm.kymograph.kgph_pvs_processed(:,fwhmloc.rem),t_axis(fwhmloc.rem))
xlim(kymograph_fig.ax,[min(t_axis(fwhmloc.rem)) max(t_axis(fwhmloc.rem))])


%%
scope = fwhmloc.nrem;
kymograph_fig = make_fig('Kymograph figure');
kymograph_fig.resolution = session.img_param.pixel2um;
kymograph_fig.update_figsize([15,3])
kymograph_fig.monitor_xyinch = [21, 5];
kymograph_fig.reset_axis
kymograph_fig.loc.x = t_axis;

kymograph_fig.plot_kymograph(session.pax_fwhm.kymograph.kgph_pvs_processed(:,scope),t_axis(scope))
xlim(kymograph_fig.ax,[min(t_axis(scope)) max(t_axis(scope))])







