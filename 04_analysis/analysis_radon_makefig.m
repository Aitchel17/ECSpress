function analysis_radon_makefig(radon_analysis, t_axis, save_dir)
% ANALYSIS_RADON_MAKEFIG Generates and saves figures for Radon analysis.
%   radon_analysis: analysis_radon object containing results.
%   t_axis: Time axis vector.
%   save_dir: Directory where figures will be saved (SVGs).

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

clee = color_lee;

% 1. Diameter Angle-Time Kymographs (Raw & Normalized)
radon_kymo_fig = make_fig("Diameter_Kymographs", 'normal');
radon_kymo_fig.update_figsize([18 3]); % 6:1 ratio for 2 subplots (each 3:1)

% Subplot 1: Raw Diameter
ax1 = subplot(1, 2, 1, 'Parent', radon_kymo_fig.fig);
radon_kymo_fig.ax = ax1;
radon_kymo_fig.plot_kymograph(radon_analysis.radon_result.diameter, t_axis, clee.gradient.inferno);
title(ax1, 'Raw Diameter');
radon_kymo_fig.put_xaxistitle('Time (s)');
radon_kymo_fig.put_yaxistitle('Angle (deg)');

% Subplot 2: Normalized Diameter
ax2 = subplot(1, 2, 2, 'Parent', radon_kymo_fig.fig);
radon_kymo_fig.ax = ax2;
radon_kymo_fig.plot_kymograph(radon_analysis.radon_result.normalized_diameterchange, t_axis, clee.gradient.inferno);
title(ax2, 'Normalized Diameter');
radon_kymo_fig.put_xaxistitle('Time (s)');
radon_kymo_fig.put_yaxistitle('Angle (deg)');

radon_kymo_fig.save2svg(save_dir);

% 2. Angle Variance (Raw & Normalized)
radon_var_fig = make_fig("Angle_Variance_Combined", 'normal');
radon_var_fig.update_figsize([10 4]);

% Subplot 1: Raw Variance
ax3 = subplot(1, 2, 1, 'Parent', radon_var_fig.fig);
radon_var_fig.ax = ax3;
radon_var_fig.plot_line(radon_analysis.radon_result.var_diameter, 'k');
title(ax3, 'Angle Variance');
radon_var_fig.put_xaxistitle('Angle (deg)');
radon_var_fig.put_yaxistitle('Variance');

% Subplot 2: Normalized Variance
ax4 = subplot(1, 2, 2, 'Parent', radon_var_fig.fig);
radon_var_fig.ax = ax4;
radon_var_fig.plot_line(radon_analysis.radon_result.var_normdiameter, 'k');
title(ax4, 'Normalized Angle Variance');
radon_var_fig.put_xaxistitle('Angle (deg)');
radon_var_fig.put_yaxistitle('Variance');

radon_var_fig.save2svg(save_dir);

% 3. Median Diameter
radonmedian_diameter_fig = make_fig("Median_diameter","normal");
radonmedian_diameter_fig.plot_line(radon_analysis.radon_result.median_diameter, 'k');
radonmedian_diameter_fig.put_xaxistitle('Angle (deg)');
radonmedian_diameter_fig.put_yaxistitle('Median Diameter (um)');
radonmedian_diameter_fig.save2svg(save_dir);

% 5. Event Images (Dilate/Constrict)
if isfield(radon_analysis.radon_result, 'events')
    for i = 1:length(radon_analysis.radon_result.events)
        event = radon_analysis.radon_result.events(i);

        % Determine if dual-channel (RGB plotting) or single channel analysis
        if size(event.median_projected, 3) == 2
            % Dual Channel: Create 1x3 Subplot (Red, Green, Merged)
            fig_title = [event.name, '_RGB_Composite'];
            radon_fig = make_fig(fig_title, 'normal');
            radon_fig.update_figsize([12 4]); % Wider figure

            img_ch1 = double(event.median_projected(:,:,1));
            img_ch2 = double(event.median_projected(:,:,2));

            % Generate RGB Image (R=Ch1, G=Ch2, B=0)
            rgb_img = plot_make_rgb(img_ch1, img_ch2);
            img_ch1_norm = rgb_img(:,:,1);
            img_ch2_norm = rgb_img(:,:,2);

            % 1. Red Channel (displayed as grayscale)
            ax1 = subplot(1, 3, 1, 'Parent', radon_fig.fig);
            radon_fig.ax = ax1;
            set(ax1, 'FontName', 'Arial', 'FontSize', 12); % Enforce Font
            imagesc(ax1, img_ch1_norm);
            axis(ax1, 'image', 'off');
            colormap(ax1, 'gray');
            title(ax1, 'Ch1 (Red)');

            center = size(img_ch1)/2 + 0.5;
            radius = min(size(img_ch1))/2 * 0.9;
            plot_radon_projection_lines(ax1, center, radius, 0:30:150, 'r');

            % 2. Green Channel (displayed as grayscale)
            ax2 = subplot(1, 3, 2, 'Parent', radon_fig.fig);
            radon_fig.ax = ax2;
            set(ax2, 'FontName', 'Arial', 'FontSize', 12); % Enforce Font
            imagesc(ax2, img_ch2_norm);
            axis(ax2, 'image', 'off');
            colormap(ax2, 'gray');
            title(ax2, 'Ch2 (Green)');
            plot_radon_projection_lines(ax2, center, radius, 0:30:150, 'r');

            % 3. Merged (RGB)
            ax3 = subplot(1, 3, 3, 'Parent', radon_fig.fig);
            radon_fig.ax = ax3;
            set(ax3, 'FontName', 'Arial', 'FontSize', 12); % Enforce Font
            imagesc(ax3, rgb_img);
            axis(ax3, 'image', 'off');
            title(ax3, 'Merged');
            plot_radon_projection_lines(ax3, center, radius, 0:30:150, 'w');

            radon_fig.save2svg(save_dir);

        else
            % Fallback: Single Channel or >2 channels generic loop
            for ch = 1:size(event.median_projected, 3)
                % Construct Figure Name
                fig_title = sprintf('%s_ch%d', event.name, ch);

                % Extract Image
                img = event.median_projected(:,:,ch);

                radon_fig = make_fig(fig_title, 'normal');
                radon_fig.showimg(img);

                center = size(img)/2 + 0.5;
                radius = min(size(img))/2 * 0.9;
                plot_radon_projection_lines(radon_fig.ax, center, radius, 0:30:150, 'r');

                radon_fig.save2svg(save_dir);
            end
        end

    end

end

disp(['Saved Radon figures to ', save_dir]);

end
