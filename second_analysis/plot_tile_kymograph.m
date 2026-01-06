function [fig] = plot_tile_kymograph(linedata_struct,opt)
%PLOT_TILE_KYMOGRAPH Summary of this function goes here
%   Detailed explanation goes here
numPlots = size(linedata_struct,2);
fig = figure('Name', 'Vessel Kymographs', 'Units', 'pixels', 'Position', [100, 100, 300, numPlots * 200]);

t = tiledlayout(round(numPlots/2), 4);

for idx = 1:numPlots
    curr_sessionid = linedata_struct(idx).sessionid;
    disp(curr_sessionid)
    kgph_bv = linedata_struct(idx).fwhmline.pax.kymograph.kgph_normbv;
    kgph_csf = linedata_struct(idx).fwhmline.pax.kymograph.kgph_normcsf;

    % bv
    ax1 = nexttile;
    imagesc(ax1, kgph_bv);
    axis(ax1, 'tight');
    axis(ax1, 'off');
    ax1.Position(3:4) = [1 0.4];  % width, height 비율 고정 예시
    title(ax1, [char(linedata_struct(idx).infodict("vesselid")),curr_sessionid ' - bv']);
    if nargin == 2
        if strcmp(opt,'csf')
            % just skip
        else
            hold on
            line_y = linedata_struct(idx).fwhmline.pax.idx.bv_lowerboundary;
            line_x = linspace(1,size(line_y,2),size(line_y,2));
            plot(ax1, line_x, line_y, 'Color', [1 0 1], 'LineWidth', 1);
            line_y = linedata_struct(idx).fwhmline.pax.idx.bv_upperboundary;
            plot(ax1, line_x, line_y, 'Color', [1 0 1], 'LineWidth', 1);
        end
        hold off
    end

    % csf
    ax2 = nexttile;
    imagesc(ax2, kgph_csf);
    axis(ax2, 'tight'); axis(ax2, 'off');
    ax2.Position(3:4) = [1 0.4];  % width, height 비율 고정 예시
    %title(ax2, [curr_field ' - csf'], 'Interpreter', 'none');
    if nargin ==2
        if strcmp(opt,'bv')
        else
                hold on
                line_y = linedata_struct(idx).fwhmline.pax.idx.pvs_lowerboundary;
                line_x = linspace(1,size(line_y,2),size(line_y,2));
                plot(ax2, line_x, line_y, 'Color', [0 1 0], 'LineWidth', 1);
                line_y = linedata_struct(idx).fwhmline.pax.idx.pvs_upperboundary;
                plot(ax2, line_x, line_y, 'Color', [0 1 0], 'LineWidth', 1);
        end
        hold off
    end

    % line_y = linedata_struct.(curr_field).fwhmline.pax.idx.bv_lowerboundary;
    % plot(ax2, line_x, line_y, 'Color', [1 0 1], 'LineWidth', 1);
    % 
    % line_y = linedata_struct.(curr_field).fwhmline.pax.idx.bv_upperboundary;
    % plot(ax2, line_x, line_y, 'Color', [1 0 1], 'LineWidth', 1);

end

t.TileSpacing = 'compact';
t.Padding = 'compact';
end

