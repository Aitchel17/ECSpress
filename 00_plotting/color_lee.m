classdef color_lee
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        gradient = struct(red = [],...
            green = [],...
            blue = [],...
            magenta = [],...
            cyan = [],...
            yellow = [],...
            gray = [],...
            inferno = [])
    end

    properties (Constant)
        clist = struct(...
            'white', [1,1,1],...
            'black', [0,0,0],...
            'cyan', [0,1,1],...
            'magenta', [1,0,1],...
            'yellow', [1,1,0],...
            'red', [1,0,0],...
            'green', [0,1,0],...
            'blue', [0,0,1],...
            'orange', [1,0.5,0],...
            'lightgreen', [0.4,0.8,0.4],...
            'darkgreen', [0,0.5,0],...
            'coloredorange', [1,0.7,0],...
            'coloredgreen', [0.6,1,0],...
            'purple', [0.5, 0, 0.5],...
            'teal', [0, 0.5, 0.5],...
            'gray', [0.5, 0.5, 0.5],...
            'lightgray', [0.8, 0.8, 0.8],...
            'darkgray', [0.2, 0.2, 0.2],...
            'brown', [0.6, 0.3, 0.1],...
            'pink', [1, 0.75, 0.8],...
            'gold', [1, 0.84, 0],...
            'navy', [0, 0, 0.5],...
            'maroon', [0.5, 0, 0]);
    end

    methods
        function obj = color_lee()
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.gradient.red =    [linspace(0,1,256)', zeros(256,1), zeros(256,1)];
            obj.gradient.green =    [zeros(256,1), linspace(0,1,256)', zeros(256,1)];
            obj.gradient.blue =    [zeros(256,1), zeros(256,1), linspace(0,1,256)'];
            obj.gradient.magenta =    [linspace(0,1,256)', zeros(256,1), linspace(0,1,256)'];
            obj.gradient.cyan =    [zeros(256,1),linspace(0,1,256)', linspace(0,1,256)'];
            obj.gradient.yellow =    [linspace(0,1,256)', linspace(0,1,256)',zeros(256,1)];
            obj.gradient.gray =    [linspace(0,1,256)', linspace(0,1,256)', linspace(0,1,256)'];

            % Inferno approximation
            x = [0, 0.15, 0.3, 0.5, 0.75, 1];
            r = [0, 0.1, 0.3, 0.7, 0.95, 1];
            g = [0, 0, 0.05, 0.2, 0.6, 1];
            b = [0, 0.2, 0.4, 0.1, 0.1, 0.8];
            points = linspace(0,1,256);
            obj.gradient.inferno = [interp1(x,r,points,'pchip')', interp1(x,g,points,'pchip')', interp1(x,b,points,'pchip')'];
            obj.gradient.inferno = max(0, min(1, obj.gradient.inferno));
        end
    end
    
    methods (Static)
        function rgb = lch(L, C, H)
            % lch Generates RGB color from Lightness, Chroma, Hue using LCH(uv) space
            % Inputs:
            %   L: Lightness (0-100)
            %   C: Chroma
            %   H: Hue (0-360)
            % vectors must be same length or scalar
            
            % 1. Create Nx1x3 image matrix for lch2rgb
            % Ensure column vectors
            if ~isscalar(L), L = L(:); end
            if ~isscalar(C), C = C(:); end
            if ~isscalar(H), H = H(:); end
            
            % Expand scalars if necessary (simple expansion)
            n_max = max([numel(L), numel(C), numel(H)]);
            if isscalar(L), L = repmat(L, n_max, 1); end
            if isscalar(C), C = repmat(C, n_max, 1); end
            if isscalar(H), H = repmat(H, n_max, 1); end
            
            lch_img = cat(3, L, C, H);
            
            % 2. Convert using lch2rgb (default mode 'luv')
            rgb_img = lch2rgb(lch_img, 'luv');
            
            % 3. Reshape back to Nx3 matrix
            rgb = squeeze(rgb_img);
            if n_max == 1
                rgb = rgb'; % Ensure row vector for single color
            end
        end

        function rgb = oklab(L, C, H)
            % oklab Generates RGB color from Lightness, Chroma, Hue using OKLAB space
             % Inputs:
            %   L: Lightness (0-1 approx, but typically scaled 0-100 in lch2rgb processing?)
            %      Ref: lch2rgb documentation says "LCH inputs ... L in [0 100]". 
            %      However, OKLAB usually uses L in [0,1]. Let's check lch2rgb implementation.
            %      lch2rgb.m line 186 defines Aok matrix. Line 191 does "imappmat(inpict,Aok/100)".
            %      This implies it expects L to be 0-100 to scale it down? 
            %      Wait, typical OKLAB L is 0-1. If user inputs 0-100, and code divides by 100, that works for 0-1 range.
            %      Let's stick to 0-100 for consistency with LCH input of lch2rgb.
            
            %   C: Chroma
            %   H: Hue (0-360)
            
             % 1. Create Nx1x3 image matrix for lch2rgb
            if ~isscalar(L), L = L(:); end
            if ~isscalar(C), C = C(:); end
            if ~isscalar(H), H = H(:); end
            
            n_max = max([numel(L), numel(C), numel(H)]);
            if isscalar(L), L = repmat(L, n_max, 1); end
            if isscalar(C), C = repmat(C, n_max, 1); end
            if isscalar(H), H = repmat(H, n_max, 1); end
            
            lch_img = cat(3, L, C, H);
            
            % 2. Convert using lch2rgb
            rgb_img = lch2rgb(lch_img, 'oklab');
            
             % 3. Reshape back to Nx3 matrix
            rgb = squeeze(rgb_img);
             if n_max == 1
                rgb = rgb'; 
            end
        end
    end

