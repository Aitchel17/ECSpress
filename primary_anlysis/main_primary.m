%% Starting of analysis pipeline
% main_primary: Load mdfExtract data
    % 1: Prepare for draw ROI and ROI based analysis
    % 1.1 Generate primary_datastruct
    % 1.2 Resample to match fps, make roianalysis struct
    % 1.3 Read roi.mat if exist
    % 1.4 generate t_axis matched with resampled 5 fps

    % 2: Draw ROI
    % 2.1: PAX for FWHM based analysis
    % 2.2: Extraparenchyma - Manually draw polygon on outter PVS boundary
    % 2.3: Dilated - Manually find dilated frame and draw polygon
    % 2.4: Constricted - Manually find constricted frame and draw polygon
    % 2.5: Dip in - Draw line start from dip in and end at center of bv
    % 2.6: Dip out - Draw line start from center and end at dip out
    % 3: Put frames of both ch1 and ch2 into roi struct, this stage
    % orginally intended to use of 
    % 4: SAVE ROILIST
 
    % 5: Generate summary figure
    % 5.1: Dilated, Constricted
            % 5.1.1 BV, constricted
            % 5.1.2 BV, dilated
            % 5.1.3 PVS, constricted
            % 5.1.4 PVS, dilated
    % 5.2: Dip in, Dip out
    % 5.3: PAX, Extraparenchyma, Dilated, Constricted
    clc,clear
clee = color_lee;
% Directory setup
directories.load_dir = 'G:\tmp\01_igkltdt\hql080\250722_hql080_whiskerb\HQL080_whiskerb_250722_007';
% directories.load_dir = 'G:\tmp\forHyunseok_251006\250502_012';

directories.save_dir = fullfile(directories.load_dir,"primary_analysis");
directories.savefig_dir =  'E:\OneDrive - The Pennsylvania State University\2023_ECSpress\08_figures';
%% 1.1 Load data
primary_datastruct = primary_loaddata(directories.load_dir);
directories.savepath = fullfile(primary_datastruct.mdfextract.info.analyzefolder,'primary_analysis');

%% 1.2 fix frame rate to 5 fps
[roianalysis.preprocessed_ch2, roianalysis.outfps_ch2] = analyze_resample(primary_datastruct.stackch2,primary_datastruct.img_param.save_fps,5);
[roianalysis.preprocessed_ch1, roianalysis.outfps_ch1] = analyze_resample(primary_datastruct.stackch1,primary_datastruct.img_param.save_fps,5);

%%
roianalysis.preprocessed_ch1 = medfilt3(roianalysis.preprocessed_ch1,[1 1 5]);
roianalysis.preprocessed_ch2 = medfilt3(roianalysis.preprocessed_ch2,[1 1 5]);
%%
% 1.3 load roi.mat if exist
if isfield(primary_datastruct,'roilist')
    roilist = primary_datastruct.roilist;
else
    roilist = roi_handle(fullfile(directories.save_dir,"roilist.mat"));
end
% 1.4 taxis generation
roianalysis.taxis = linspace(primary_datastruct.img_param.imgstarttime,...
    primary_datastruct.img_param.imgendtime,size(roianalysis.preprocessed_ch1,3));
% resolution
primary_datastruct.img_param.pixel2um


%% 2.1 PAX :  Run analysis_pax_fwhm
% Primary axis (pax) - Where successful Kymograph generation to be guarented
try
roilist.addroi(roianalysis.preprocessed_ch1,'pax','line');
catch ME
    disp(ME.message)
    roilist.modifyroi(roianalysis.preprocessed_ch2,'pax');
end


%% 2.2 Extraparenchyma : 
% Extraparenchyma - For rough calculation of constricted, dilated state
try
roilist.addroi(roianalysis.preprocessed_ch2,'extraparenchyma','polygon');
catch ME
    disp(ME.message)
    roilist.modifyroi(roianalysis.preprocessed_ch2,'extraparenchyma');
end
%% 2.3 Constricted bv
try
    roilist.copyroi('extraparenchyma','constricted_bv');
catch ME
    disp(ME.message)
    roilist.modifyroi(roianalysis.preprocessed_ch1,'constricted_bv');
end

%% 2.4 Dilated bv
try
    roilist.copyroi('constricted_bv','dilated_bv');
catch ME
    disp(ME.message)
    roilist.modifyroi(roianalysis.preprocessed_ch1,'dilated_bv');
end
%% 2.5 Dip in angle : 
% Extraparenchyma - For rough calculation of constricted, dilated state
try
roilist.addroi(roianalysis.preprocessed_ch1,'dipin','line');
catch ME
    disp(ME.message)
    roilist.modifyroi(roianalysis.preprocessed_ch1,'dipin');
end
%% 2.6 Dip out angle : 
% Extraparenchyma - For rough calculation of constricted, dilated state
try
roilist.addroi(roianalysis.preprocessed_ch2,'dipout','line');
catch ME
    disp(ME.message)
    
    roilist.modifyroi(roianalysis.preprocessed_ch1,'dipout');
end

%%  2.6 radon : 
% Extraparenchyma - For rough calculation of constricted, dilated state
try
roilist.addroi(roianalysis.preprocessed_ch1,'radon','rectangle');
catch ME
    disp(ME.message)
    
    roilist.modifyroi(roianalysis.preprocessed_ch2,'radon');
end

%%

%% 3. Put image slice for Level 1 ROI list
roilist.addimgchannel(cat(4,roianalysis.preprocessed_ch1,roianalysis.preprocessed_ch2),...
    ["pax","extraparenchyma","dipin","dipout","dilated_bv","constricted_bv"]); % put image slices

%% 4. SAVE ROILIST
roilist.save2disk

%% 5.1.1 BV_constricted + dilated and constricted vessel outline
fig.roi_dilatedbv = make_fig('roi_dilated_bv');
fig.roi_dilatedbv.update_figsize([8 6])
fig.roi_dilatedbv.reset_axis
fig.roi_dilatedbv.showrois(roilist,1,["dilated_bv","constricted_bv"],["-r","-g"])
%% Adjust contrast and save data
fig.roi_dilatedbv.save2svg(directories.save_dir) % save
%% 5.1.2 BV_constricted + dilated and constricted vessel outline

fig.roi_constrictedbv = make_fig('roi_constricted_bv');
fig.roi_constrictedbv.update_figsize([8 6])
fig.roi_constrictedbv.reset_axis
fig.roi_constrictedbv.showrois(roilist,1,["constricted_bv","dilated_bv"],["-g","-r"])
%% Adjust contrast and save data
fig.roi_constrictedbv.save2svg(directories.save_dir) % save
%% 5.1.3 PVS_dilated image + dilated and constricted vessel outline
fig.roi_dilatedpvs = make_fig('roi_dilated_pvs');
fig.roi_dilatedpvs.update_figsize([8 6])
fig.roi_dilatedpvs.reset_axis
fig.roi_dilatedpvs.showrois(roilist,2,["dilated_bv","constricted_bv"],["-r","-g"])
%% Adjust contrast and save data
fig.roi_dilatedpvs.save2svg(directories.save_dir)
%% 5.1.4 PVS_constricted image + dilated and constricted vessel outline
fig.roi_constrictedpvs = make_fig('roi_constricted_pvs');
fig.roi_constrictedpvs.update_figsize([8 6])
fig.roi_constrictedpvs.reset_axis
fig.roi_constrictedpvs.showrois(roilist,2,["constricted_bv","dilated_bv"],["-g","-r"])
%%
fig.roi_constrictedpvs.save2svg(directories.save_dir)
%% 5.2 Dip in, Dip out
fig.roi_dipinout = make_fig('roi_dipinout');
fig.roi_dipinout.update_figsize([8 6])
fig.roi_dipinout.reset_axis
fig.roi_dipinout.showrois(roilist,1,["dipin","dipout"],["-g","-r"])
%%
fig.roi_dipinout.save2svg(directories.save_dir);


%% 5.3 PAX, Dilated, Constricted
fig.roi_pax = make_fig('roi_pax');
fig.roi_pax.update_figsize([8 6])
fig.roi_pax.bring_fig
fig.roi_pax.reset_axis
fig.roi_pax.showrois(roilist,3,["pax","dipin","dipout","constricted_bv","dilated_bv"],["-c","-y","-y","y","y"])
%%
fig.roi_pax.save2svg(directories.save_dir);
%%
fig.roi_pax = make_fig('roi_pax');
fig.roi_pax.update_figsize([8 6])
fig.roi_pax.bring_fig
fig.roi_pax.reset_axis
fig.roi_pax.showrois(roilist,2,"pax",'-c')
