% folder directory structure
    % Mouse ID
        % Imaging date_mouseID_expClass* (sleep, whisker, whiskerball, ball)
            % mouseID_expClassdate_fileid
clc,clear
folderpath ='G:\tmp\hql072\**';

dir_list = struct2table(dir(folderpath));

line_dirloc = matches(dir_list.name,'line_fwhms.mat');
line_list = dir_list(line_dirloc,:);

info_dirloc = endsWith(dir_list.name, '_info.txt');
info_list = dir_list(info_dirloc,:);

analog_dirloc = matches(dir_list.name, 'analog.mat');
analog_list = dir_list(analog_dirloc,:);



% path
linedata_struct = struct();
for i = 1:height(line_list)
    linedata_path = fullfile(line_list.folder{i},line_list.name{i});
    
    info_listloc = matches(info_list.folder, fileparts(line_list.folder{i}));
    info_dir = info_list(info_listloc,:);
    info_path = fullfile(info_dir.folder{1},info_dir.name{1});
    disp(info_path)
    analog_listloc = matches(analog_list.folder,line_list.folder{i});
    analog_dir = analog_list(analog_listloc,:);
    analog_path = fullfile(analog_dir.folder{1},analog_dir.name{1});
    
    
    info_table = readtable(info_path);
    mdforigin = info_table(matches(info_table.Field,'mdfName'),:).Value{1};
    fieldname = mdforigin(1:end-4);
    disp(fieldname)

    load(analog_path)
    load(linedata_path)
    

    linedata_struct.(fieldname).analog = load(analog_path).primary_analog;
    linedata_struct.(fieldname).fwhmline = load(linedata_path).line_fwhms;
    linedata_struct.(fieldname).infotable = info_table;

end
%%
sep_thickness = 2; % 구분선 두께
sep_color = 1;     % 흰색: 1, 검정색: 0

fieldnames_top = fieldnames(linedata_struct);
numPlots = numel(fieldnames_top);
tileCols = ceil(sqrt(numPlots));
tileRows = ceil(numPlots / tileCols);

figure('Name', 'Vessel Kymographs');
t = tiledlayout(tileRows, tileCols, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:numPlots
    nexttile;

    curr_field = fieldnames_top{i};
    kgph_bv = linedata_struct.(curr_field).fwhmline.pax.kymograph.kgph_normbv;
    kgph_csf = linedata_struct.(curr_field).fwhmline.pax.kymograph.kgph_normcsf;

    % 두 이미지 사이에 구분선 삽입
    separator = sep_color * ones(sep_thickness, size(kgph_bv, 2));
    combined_img = [kgph_bv; separator; kgph_csf];

    imagesc(combined_img);
    axis tight;
    axis off;
    title(curr_field, 'Interpreter', 'none');
end
%%
fieldnames_top = fieldnames(linedata_struct);
numPlots = numel(fieldnames_top);

figure('Name', 'Vessel Kymographs', 'Units', 'pixels', 'Position', [100, 100, 300, numPlots * 200]);

t = tiledlayout(round(numPlots/2), 4);

for i = 1:numPlots
    curr_field = fieldnames_top{i};

    kgph_bv = linedata_struct.(curr_field).fwhmline.pax.kymograph.kgph_normbv;
    kgph_csf = linedata_struct.(curr_field).fwhmline.pax.kymograph.kgph_normcsf;

    % bv
    ax1 = nexttile;
    imagesc(ax1, kgph_bv);
    axis(ax1, 'tight'); axis(ax1, 'off');
    ax1.Position(3:4) = [1 0.4];  % width, height 비율 고정 예시
    title(ax1, [curr_field ' - bv'], 'Interpreter', 'none');
    hold on
    line_y = linedata_struct.(curr_field).fwhmline.pax.idx.bv_lowerboundary;
    line_x = linspace(1,size(line_y,2),size(line_y,2));
    plot(ax1, line_x, line_y, 'Color', [1 0 1], 'LineWidth', 1);
    line_y = linedata_struct.(curr_field).fwhmline.pax.idx.bv_upperboundary;
    plot(ax1, line_x, line_y, 'Color', [1 0 1], 'LineWidth', 1);
    % csf
    ax2 = nexttile;
    imagesc(ax2, kgph_csf);
    axis(ax2, 'tight'); axis(ax2, 'off');
    ax2.Position(3:4) = [1 0.4];  % width, height 비율 고정 예시
    %title(ax2, [curr_field ' - csf'], 'Interpreter', 'none');
    hold on
    line_y = linedata_struct.(curr_field).fwhmline.pax.idx.pvs_lowerboundary;
    line_x = linspace(1,size(line_y,2),size(line_y,2));
    plot(ax2, line_x, line_y, 'Color', [0 1 0], 'LineWidth', 1);
    line_y = linedata_struct.(curr_field).fwhmline.pax.idx.pvs_upperboundary;
    plot(ax2, line_x, line_y, 'Color', [0 1 0], 'LineWidth', 1);
    % line_y = linedata_struct.(curr_field).fwhmline.pax.idx.bv_lowerboundary;
    % plot(ax2, line_x, line_y, 'Color', [1 0 1], 'LineWidth', 1);
    % 
    % line_y = linedata_struct.(curr_field).fwhmline.pax.idx.bv_upperboundary;
    % plot(ax2, line_x, line_y, 'Color', [1 0 1], 'LineWidth', 1);

end

t.TileSpacing = 'compact';
t.Padding = 'compact';
%%
fieldnames_top = fieldnames(linedata_struct);
numPlots = numel(fieldnames_top);

figure('Name', 'Vessel Kymographs', 'Units', 'pixels', 'Position', [100, 100, 300, numPlots * 200]);

t = tiledlayout(round(numPlots/2), 4);

for i = 1:numPlots
    curr_field = fieldnames_top{i};

    kgph_bv = linedata_struct.(curr_field).fwhmline.pax.kymograph.kgph_normbv;
    kgph_csf = linedata_struct.(curr_field).fwhmline.pax.kymograph.kgph_normcsf;

    % bv
    ax1 = nexttile;
    imagesc(ax1, kgph_bv);
    axis(ax1, 'tight'); axis(ax1, 'off');
    ax1.Position(3:4) = [1 0.4];  % width, height 비율 고정 예시
    title(ax1, [curr_field ' - bv'], 'Interpreter', 'none');

    % csf
    ax2 = nexttile;
    imagesc(ax2, kgph_csf);
    axis(ax2, 'tight'); axis(ax2, 'off');
    ax2.Position(3:4) = [1 0.4];  % width, height 비율 고정 예시
    %title(ax2, [curr_field ' - csf'], 'Interpreter', 'none');


end

t.TileSpacing = 'compact';
t.Padding = 'compact';


%% Draw horizontal line on top
fig = figure;
ax = axes('Parent', fig);
line_y = linedata_struct.(curr_field).fwhmline.pax.idx.pvs_lowerboundary;
line_x = linspace(1,size(line_y,2),size(line_y,2));
plot(ax1, line_x, line_y, 'Color', [0 1 0], 'LineWidth', 1);
line_y = linedata_struct.(curr_field).fwhmline.pax.idx.pvs_upperboundary;
plot(ax1, line_x, line_y, 'Color', [0 1 0], 'LineWidth', 1);
%%




























































 