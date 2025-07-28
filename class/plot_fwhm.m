classdef plot_fwhm
    %PLOT_FWHM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = plot_fwhm(inputArg1,inputArg2)
            %PLOT_FWHM Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %% kymograph with overlay
            fig = figure;
            ax = axes('Parent', fig);
            
            % Simulated kymograph data
            kymograph = pax_fwhm.kymograph.kgph_csf;
            kymograph = (kymograph-min(kymograph,[],1))./max(kymograph,[],1);
            
            % Display kymograph using imagesc with colormap
            imagesc(ax, kymograph);
            colormap(ax, parula);  % Or your inferno colormap
            hold(ax, 'on');
            
            % Draw horizontal line on top
            line_y = pax_fwhm.idx.pvs_lowerboundary;
            line_x = linspace(1,size(line_y,2),size(line_y,2));
            plot(ax, line_x, line_y, 'Color', [0 1 0], 'LineWidth', 1);
            line_y = pax_fwhm.idx.pvs_upperboundary;
            plot(ax, line_x, line_y, 'Color', [0 1 0], 'LineWidth', 1);
            line_y = pax_fwhm.idx.bv_lowerboundary;
            plot(ax, line_x, line_y, 'Color', [1 0 1], 'LineWidth', 1);
            line_y = pax_fwhm.idx.bv_upperboundary;
            plot(ax, line_x, line_y, 'Color', [1 0 1], 'LineWidth', 1);
        end
    end
end

