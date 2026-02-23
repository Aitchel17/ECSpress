clc, clear
clee = color_lee();
exp_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis\00_igkl';
table_path = fullfile(exp_path, "transition.mat");
load_struct = load(table_path);
transition_table = load_struct.save_content;
%%


cred = clee.lch(35,150,15);

cgreen = clee.lch(35,60,140);
cblue = clee.lch(35,100,260);
%%

sess_ave = transition_table.awakenorm_data.Date_ave;
mouse_ave = transition_table.awakenorm_data.MouseID_ave;



%%
logic.state = sess_ave.state_name == "nr_trans";
logic.dtype = sess_ave.DataType == "thickness_bv";
combined_logic = all(table2array(struct2table(logic)),2);
nr_trans.bvthickness = sess_ave(combined_logic,:);
%%
logic = struct();
logic.state = mouse_ave.state_name == "nr_trans";
logic.dtype = mouse_ave.DataType == "thickness_bv";
combined_logic = all(table2array(struct2table(logic)),2);
nr_trans.mbvthickness = mouse_ave(combined_logic,:);


%%
cla
cfaintred = clee.lch(100,130,10);
cred = clee.lch(35,150,15);
plot(cell2mat(nr_trans.bvthickness.data)',"Color",cfaintred);
hold on
plot(cell2mat(nr_trans.mbvthickness.data),"Color",cred)

xlim([0 150])
ylim([0.8 1.5])

%%
figure()

