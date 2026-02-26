clc, clear
addpath(genpath('g:\03_program\01_ecspress\00_plotting'));
clee = color_lee();
exp_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis\00_igkl';
table_path = fullfile(exp_path, "transition.mat");
load_struct = load(table_path);
transition_table = load_struct.save_content;

%% Shared definitions
transitions  = ["an_trans", "na_trans", "nr_trans", "ra_trans"];
data_types = ["thickness_bv", "thickness_totalpvs", "thickness_eps"];
%
figconfig = struct();

figconfig.bv.faint  = clee.lch(80, 130,  10);
figconfig.bv.bold   = clee.lch(45, 150,  15);

figconfig.pvs.faint = clee.lch(80, 130, 110);
figconfig.pvs.bold  = clee.lch(45, 150, 115);

figconfig.eps.faint = clee.lch(80, 130,  70);
figconfig.eps.bold  = clee.lch(45, 120,  75);
%
vessel_names = ["BV", "PVS", "EPS"];
vessel_colors = {figconfig.bv, figconfig.pvs, figconfig.eps};

label_map = struct( ...
    'an_trans', struct('pre',"Awake",'post',"NREM"), ...
    'na_trans', struct('pre',"NREM", 'post',"Awake"), ...
    'nr_trans', struct('pre',"NREM", 'post',"REM"), ...
    'ra_trans', struct('pre',"REM",  'post',"Awake") ...
);

%%
norm_fields = ["abs_numeric", "awakenorm_numeric", "awakesubtract_numeric"];
transitions  = ["an_trans", "na_trans", "nr_trans", "ra_trans"];
transition_labels = ["Awake2NREM", "NREM2Awake", "NREM2REM", "REM2Awake"];
data_types = ["thickness_bv", "thickness_totalpvs", "thickness_eps"];
plotdata_dtypes = ["bv","totalpvs","eps"];
%%
sess_ave  = transition_table.(norm_fields(2)).Date_ave;
mouse_ave = transition_table.(norm_fields(2)).MouseID_ave;
plotdata  = build_plotdata(sess_ave, mouse_ave, transitions, data_types);
%%

%%
plot_data.meanchanges = [];
plot_data.stdchanges = [];
plot_data.xlabel = [];
for tname_idx = 1:numel(transitions)
    for dname_idx = 1:numel(plotdata_dtypes)
        tmp.dname = plotdata_dtypes(dname_idx);
        %tname_idx = 1;
        tmp.tname = transitions(tname_idx);
        tmp.tlabel = transition_labels(tname_idx);
        tmp.labelname = strcat(tmp.tlabel," ",tmp.dname);
        tmp.changes = plotdata.(tmp.tname).(tmp.dname).sess.post_median-plotdata.(tmp.tname).(tmp.dname).sess.pre_median;
        plot_data.meanchanges = [plot_data.meanchanges mean(tmp.changes)];
        plot_data.stdchanges = [plot_data.stdchanges std(tmp.changes,0,"all")];
        plot_data.xlabel = [plot_data.xlabel tmp.labelname];
    end
end
%%
plot_data.meanchanges = [];
plot_data.stdchanges = [];
plot_data.xlabel = [];
for dname_idx = 1:numel(plotdata_dtypes)
    for tname_idx = 1:numel(transitions)
        tmp.dname = plotdata_dtypes(dname_idx);
        %tname_idx = 1;
        tmp.tname = transitions(tname_idx);
        tmp.tlabel = transition_labels(tname_idx);
        tmp.labelname = strcat(tmp.tlabel," ",tmp.dname);
        tmp.changes = plotdata.(tmp.tname).(tmp.dname).sess.post_median-plotdata.(tmp.tname).(tmp.dname).sess.pre_median;
        plot_data.meanchanges = [plot_data.meanchanges mean(tmp.changes)];
        plot_data.stdchanges = [plot_data.stdchanges std(tmp.changes,0,"all")];
        plot_data.xlabel = [plot_data.xlabel tmp.labelname];
    end
end


%%
figure()
cla, clc
bar(plot_data.xlabel,plot_data.meanchanges)
hold on
errorbar(plot_data.xlabel,plot_data.meanchanges,plot_data.stdchanges)
%%


function plotdata = build_plotdata(sess_tbl, mouse_tbl, transitions, datatypes)
% BUILD_PLOTDATA  Index all data by (transition, datatype).
%   plotdata.(trans).(dtype).sess  = filtered session-level table rows
%   plotdata.(trans).(dtype).mouse = filtered mouse-level table rows
    for ti = 1:numel(transitions)
        tr = transitions(ti);
        for di = 1:numel(datatypes)
            dt = datatypes(di);
            % field name: replace non-identifier chars
            dtfield = strrep(dt, "thickness_", "");
            plotdata.(tr).(dtfield).sess  = ...
                sess_tbl(sess_tbl.state_name == tr & sess_tbl.DataType == dt, :);
            plotdata.(tr).(dtfield).mouse = ...
                mouse_tbl(mouse_tbl.state_name == tr & mouse_tbl.DataType == dt, :);
        end
    end
end
