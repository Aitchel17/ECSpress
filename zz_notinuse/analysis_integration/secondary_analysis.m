% folder directory structure
    % Mouse ID
        % Imaging date_mouseID_expClass* (sleep, whisker, whiskerball, ball)
            % mouseID_expClassdate_fileid
clee = color_lee;
%%
% To-do:
% 1. Successful kymograph selection (for a single animal)
% 2. Construct heatmap data structure.
%    Required for heatmap generation:
%    2.1. Datasets: vascular change vs. total PVS change,
%         vascular change vs. unidirectional PVS change,
%         vascular change vs. other-direction PVS change
%    2.2. Scale and magnification differ between sessions...

% 3. Merge heatmap data
%    3.1. Data required for merging: origin points — mode of vascular change and mode of PVS change
%    3.2. Array of PVS mode values corresponding to each vascular change
%    3.3. Origin alignment (mode-based)
%    3.4. 95% confidence interval values
%    3.5. Patch-based CI visualization and linear piecewise fitting

% 5. Fit slopes (angles) of x and y axes — first segment angle and second segment angle
%%
clc, clear
%%
foldername = '01_igkltdt';
mouseid = 'hql080';
path = fullfile('G:\tmp\', foldername, mouseid, '**');
%%
savepath= struct();
savepath.loadsplit = strsplit(path,filesep);
savepath.primarystruct = fullfile(savepath.loadsplit{1:end-1},'primarystruct.mat');
savepath.secondarystruct = fullfile(savepath.loadsplit{1:end-1},'paxfwhm_struct.mat');
%%
primarystruct = primary_integration(string(path));
%% 1. kymograph selection
tmp.fig = plot_tile_kymograph(primarystruct,'all');
ui_tile_kymograph_selection(tmp.fig,{primarystruct.sessionid})
%%
tmp.goodsessionlist = getappdata(tmp.fig, 'clicked_list');
secondary_struct = repmat(struct(),1,length(tmp.goodsessionlist)); 
thickness_fields = {'bv','uppvs', 'downpvs', 'totalpvs'};
heat_types = {'uppvs', 'downpvs', 'totalpvs'};
for idx = 1:length(tmp.goodsessionlist)
    % 세션 불러오기
    secondary_struct(idx).session_id = tmp.goodsessionlist(idx);
    [~,tmp.sessionloci] = ismember(tmp.goodsessionlist,{primarystruct.sessionid});
    primary_session = primarystruct(tmp.sessionloci(idx));
    secondary_session = struct(); 
    % 이름 정의
    tmp.name_parts = strsplit(primary_session.infodict("mdfName"), '.');
    tmp.name = tmp.name_parts{1};

    % 후처리
    tmp.scale_parts = strsplit(primary_session.infodict("objpix"));
    tmp.scale = str2double(tmp.scale_parts{1});
   
    tmp.thickness = kymograph_afterprocessing(primary_session);
    for tidx = 1:4
    secondary_session.thickness(tidx) = struct( ...
        'type', thickness_fields{tidx}, ...
        'thickness', tmp.thickness.([thickness_fields{tidx} '_thickness']), ...
        'changes', tmp.thickness.([thickness_fields{tidx} '_changes']));
    end

    % heatdata 초기화
    secondary_session.heatdata = repmat(struct(), 1, 3);
    
    for hidx = 1:3
        secondary_session.heatdata(hidx).type = heat_types{hidx};
        
        % heatmap 임시 데이터
        tmp.heatdata = xy2heatmap(...
            secondary_session.thickness(1).changes,...
            secondary_session.thickness(hidx+1).changes,...
            tmp.scale);
        secondary_session.heatdata(hidx).xy_counts = tmp.heatdata.xy_counts;
        secondary_session.heatdata(hidx).x_centers = tmp.heatdata.x_centers;
        secondary_session.heatdata(hidx).y_centers = tmp.heatdata.y_centers;

        % 히트맵 후처리 결과를 곧바로 병합하여 heatdata에 저장
        tmp_heatpost = heatmap_postprocessing(tmp.heatdata);
        secondary_session.heatdata(hidx).xy_counts_aligned = tmp_heatpost.xy_counts_clean;
        secondary_session.heatdata(hidx).x_centers_aligned = tmp_heatpost.x_baseceneters;
        secondary_session.heatdata(hidx).y_centers_aligned = tmp_heatpost.y_baseceneters;

        secondary_session.heatdata(hidx).modespvs = tmp_heatpost.modepvs;
        secondary_session.heatdata(hidx).logxycount = tmp.heatdata.log_xycounts;
    end
    if length(secondary_session.heatdata(1).y_centers_aligned) <length(secondary_session.heatdata(2).y_centers_aligned)
        disp('uppvs and downpvs location flip')
        secondary_session.heatdata = secondary_session.heatdata([2,1,3]);
    end
        secondary_struct(idx).thickness = secondary_session.thickness;
        secondary_struct(idx).heatdata = secondary_session.heatdata;
end

%%
save(savepath.primarystruct,'primarystruct')
secondary_paxfwhm = secondary_struct;
save(savepath.secondarystruct,'secondary_paxfwhm')

 