classdef make_fig < handle
    %MAKE_FIG Summary of this class goes here
    %   Detailed explanation goes here

    properties
        save_path = 'E:\OneDrive - The Pennsylvania State University\2023_ECSpress\08_figures';
        monitor_xyinch = [27 2];
        xy_sizeinch = [5 2];
        fig_color = 'w';
        axis_color = 'k';
        fontsize = 12;
        fontname = 'Arial';
        resolution = 1; % Default y axis unit: pixels
        fps = 1; % Default x axis unit:frame
        blackbackground = 0;
        fig
        axis_type
        ax
        loc = struct('x',[],'y',[]);
    end

    methods
        function obj = make_fig(fig_name,axis_opt)
            arguments
                fig_name char
                axis_opt char {mustBeMember(axis_opt, {'normal','polar'})} = 'normal'
            end
            obj.fig = figure("Name",fig_name);
            set(obj.fig,'Units','inches',...
                "Position",[obj.monitor_xyinch(1) obj.monitor_xyinch(2) obj.xy_sizeinch(1) obj.xy_sizeinch(2)])
            obj.axis_type = axis_opt;
            if strcmp(axis_opt, 'normal')
                obj.initialize_axis;
            elseif strcmp(axis_opt, 'polar')
                obj.ax = polaraxes(obj.fig);
            end
        end

        function plot_polar(obj,theta,angular_position,color,marker)
            % angular_position = [angular_position,angular_position(1)];
            % theta = linspace(0,2*pi,length(angular_position));
            pplot = polarplot(obj.ax,theta,angular_position);
            pplot.Color = color;
            pplot.Marker = marker;
        end

        function polar_init(obj,angular_position,color,marker)
            cla(obj.ax)
        end

        function plot_kymograph(obj,kymograph_data,taxis,color)
            arguments
                obj
                kymograph_data
                taxis = []
                color = 'gray'
            end
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            numel_x = size(kymograph_data,2);
            numel_y = size(kymograph_data,1);
            if ~isempty(taxis)
                obj.loc.x = taxis;
            else
                obj.loc.x = linspace(0,(numel_x-1)/obj.fps,numel_x); % start from 0
            end
            obj.loc.y = linspace(0,(numel_y-1)*obj.resolution,numel_y); % start from 0
            hold(obj.ax,"on")
            imagesc(obj.ax,obj.loc.x, obj.loc.y, kymograph_data);
            xlim(obj.ax,[0 max(obj.loc.x)]);
            ylim(obj.ax,[0 max(obj.loc.y)]);
            colormap(color)
        end



        function plot_line(obj,line_data,color,marker, linestyle, linewidth)
            arguments
                obj
                line_data
                color
                marker = 'none'
                linestyle = '-'
                linewidth = 1
            end
            numel_x = length(line_data);
            if isempty(obj.loc.x)
                obj.loc.x = linspace(0,(numel_x-1)/obj.fps,numel_x);
            else
                if numel_x == numel(obj.loc.x)
                    disp('using preregistered x,y axis location information')
                else
                    warning('The line_data length is %d while preregistered axis length is %d, reset ',...
                        numel_x, numel(obj.loc.x))
                    obj.loc.x = linspace(0,(numel_x-1)/obj.fps,numel_x);
                end
            end
            p = plot(obj.ax,obj.loc.x,line_data*obj.resolution);
            p.Color = color;
            p.Marker = marker;
            p.LineStyle =  linestyle;
            p.LineWidth = linewidth;
            xlim(obj.ax,[0 max(obj.loc.x)]);
            if isempty(obj.loc.y)
                ylim(obj.ax,[min(line_data*obj.resolution) max(line_data*obj.resolution)]);
            else
                ylim(obj.ax,[0 max(obj.loc.y)]);
            end
            obj.preset_axis
        end

        function plot_xline(obj,line_coordination,line_spec)
            if ~isempty(obj.loc.x)
                line_coordination = obj.loc.x(line_coordination);
            end
            for k = 1:length(line_coordination)
                xline(obj.ax,line_coordination-0.5,line_spec,'LineWidth',1); % x axis start position should be 0
            end
        end

        function plot_scatter(obj,x_axis,y_axis,scatter_spec)
            obj.loc.x = x_axis;
            obj.loc.y = y_axis;
            min_xy = min(min(x_axis,y_axis));
            max_xy = max(max(x_axis,y_axis));

            scatter(obj.ax,obj.loc.x,obj.loc.y,scatter_spec)
            hold on
            plot([min_xy, max_xy], [min_xy, max_xy], 'w')
            xlim(obj.ax,[min_xy max_xy]);
            ylim(obj.ax,[min_xy max_xy]);
            obj.preset_axis
        end


    end
    %%  ROI visualization related methods
    methods
        function showroi(obj, roi_handle, label, channel)
            arguments
                obj
                roi_handle
                label
                channel
            end
            i = roi_handle.findLabel(label);
            image = roi_handle.ROIs(i).ROISlice;
            vertices = roi_handle.ROIs(i).Vertices;
            roimod = roi_handle.ROIs(i).Mode;

            %%
            im_handle = imagesc(obj.ax,image(:,:,channel)); % plot image
            axis(obj.ax, 'image')
            colormap('gray')

            hold on;

            if strcmp(roimod, 'line')
                plot(obj.ax,[vertices(1, 1); vertices(2, 1)], [vertices(1, 2); vertices(2, 2)], 'r-', 'LineWidth', vertices(3, 1)); % Close the polygon
            else
                plot(obj.ax,[vertices(:, 1); vertices(1, 1)], [vertices(:, 2); vertices(1, 2)], 'r-', 'LineWidth', 2); % Close the polygon
            end
            hold off;

            % Output the original image (no modifications made to the image itself)
            imcontrast(im_handle)
        end

        function showimg(obj,image)
            im_handle = imagesc(obj.ax,image(:,:)); % plot image
            axis(obj.ax, 'image')
            colormap('gray')
            % Output the original image (no modifications made to the image itself)
            imcontrast(im_handle)
        end

        function showrois(obj, roi_handle, channel, labels,colorlist) % 250925, display all region of interest
            arguments
                obj
                roi_handle
                channel
                labels = []
                colorlist = []
            end
            if isempty(labels)
                image = zeros(size(roi_handle.ROIs(1).ROISlice));
            else
                i = roi_handle.findLabel(labels(1));
                image = roi_handle.ROIs(i).ROISlice;
            end
            %%
            idx_list = [];
            for i = 1:length(labels)
                idx = roi_handle.findLabel(labels(i));
                idx_list = [idx_list, idx];
            end
            % 251015_rgb mode added
            if channel == 3
                rgb_image = plot_make_rgb(image(:,:,1), image(:,:,2));
                im_handle = imagesc(obj.ax,rgb_image); % plot image
            else
                im_handle = imagesc(obj.ax,image(:,:,channel)); % plot image
            end

            axis(obj.ax, 'image')
            colormap('gray')
            hold on
            for i = 1:length(idx_list)
                vertices = roi_handle.ROIs(idx_list(i)).Vertices;
                roimod = roi_handle.ROIs(idx_list(i)).Mode;
                lineprofile = colorlist(i);
                if strcmp(roimod, 'line')
                    plot(obj.ax,[vertices(1, 1); vertices(2, 1)], [vertices(1, 2); vertices(2, 2)], lineprofile, 'LineWidth', vertices(3, 1)); % Close the polygon
                else
                    plot(obj.ax,[vertices(:, 1); vertices(1, 1)], [vertices(:, 2); vertices(1, 2)], lineprofile, 'LineWidth', 1); % Close the polygon
                end
            end

            % Output the original image (no modifications made to the image itself)
            if channel ~=3
                imcontrast(im_handle)
            end
        end

        function save2svg(obj,save_path)
            arguments
                obj
                save_path = [];
            end
            %%
            if isempty(save_path)
                path = fullfile(obj.save_path, obj.fig.Name);
            else
                path = fullfile(save_path, obj.fig.Name);
            end
            print(obj.fig, path, "-dsvg", "-vector");

        end

    end


    %% axis
    methods
        function initialize_axis(obj)
            if isempty(obj.ax)
                obj.ax = axes(obj.fig);
            end
            obj.preset_axis

        end

        function update_figsize(obj,size_inch)
            obj.xy_sizeinch = size_inch;
            set(obj.fig,'Units','inches',...
                "Position",[obj.monitor_xyinch(1) obj.monitor_xyinch(2) obj.xy_sizeinch(1) obj.xy_sizeinch(2)])
        end

        function update_position(obj,position_inch)
            obj.monitor_xyinch = position_inch;
            set(obj.fig,'Units','inches',...
                "Position",[obj.monitor_xyinch(1) obj.monitor_xyinch(2) obj.xy_sizeinch(1) obj.xy_sizeinch(2)])
        end

        function reset_axis(obj)
            cla(obj.ax);
            obj.loc = struct('x',[],'y',[]);
        end


        function bring_fig(obj)
            figure(obj.fig)
        end

        function hold_axis(obj,bool)
            if bool
                hold(obj.ax,'on')
            else
                hold(obj.ax,'off')
            end
        end

        function convert_background(obj,bool_blackbackground)
            obj.blackbackground = bool_blackbackground;
            if obj.blackbackground
                obj.fig_color = 'k';
                obj.axis_color = 'w';
            else
                obj.fig_color = 'w';
                obj.axis_color = 'k';
            end

            obj.preset_axis
        end
        function put_xaxistitle(obj,xaxis_title)
            obj.ax.XLabel.String = xaxis_title;
            obj.preset_axis
        end

        function put_yaxistitle(obj,yaxis_title)
            obj.ax.YLabel.String = yaxis_title;
            obj.preset_axis
        end

        function change_xylim(obj, xrange, yrange)
            arguments
                obj
                xrange = obj.loc.x
                yrange = obj.loc.y
            end
            xlim(obj.ax,xrange);
            ylim(obj.ax,yrange);
        end
    end
    methods (Access = private)
        function preset_axis(obj)

            % black background related setup
            set(obj.fig,'Color',obj.fig_color)
            set(obj.ax, 'Color',obj.fig_color, ...
                'XColor', obj.axis_color, ...
                'YColor', obj.axis_color);
            % general axis setup
            set(obj.ax,'Box','off', ...
                'TickDir','out', ...
                'LineWidth',1, ...
                'FontName',obj.fontname, ...
                'FontSize',obj.fontsize, ...
                'Layer','top')
        end
    end


end

