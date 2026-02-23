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
            % lch Converts LCH(uv) color to sRGB.
            % L: Lightness 0-100, C: Chroma, H: Hue 0-360 degrees
            % Returns Nx3 matrix of clamped sRGB values [0,1]
            
            % Ensure column vectors
            if ~isscalar(L), L = L(:); end
            if ~isscalar(C), C = C(:); end
            if ~isscalar(H), H = H(:); end
            n = max([numel(L), numel(C), numel(H)]);
            if isscalar(L), L = repmat(L, n, 1); end
            if isscalar(C), C = repmat(C, n, 1); end
            if isscalar(H), H = repmat(H, n, 1); end

            % LCH(uv) -> LUV
            Hrad = H * pi / 180;
            U = C .* cos(Hrad);
            V = C .* sin(Hrad);

            % LUV -> XYZ  (D65 white point)
            WP = [0.95047, 1.0, 1.08883];
            refd = WP(1) + 15*WP(2) + 3*WP(3);
            refU = 4*WP(1) / refd;
            refV = 9*WP(2) / refd;

            fY = (L + 16) / 116;
            ep = 216/24389; kp = 24389/27;
            Ylin = fY.^3;
            mask = Ylin < ep;
            Ylin(mask) = (116*fY(mask) - 16) / kp;

            mk = (L == 0);
            Un = U ./ (13*L + 1e-10*mk) + refU;
            Vn = V ./ (13*L + 1e-10*mk) + refV;

            Xn = -(9 * Ylin .* Un) ./ ((Un - 4) .* Vn - Un .* Vn);
            Zn = (9*Ylin - 15*Vn.*Ylin - Vn.*Xn) ./ (3*Vn);

            % XYZ -> linear sRGB  (D65 sRGB matrix)
            M = [3.2404542, -1.5371385, -0.4985314;
                -0.9692660,  1.8760108,  0.0415560;
                 0.0556434, -0.2040259,  1.0572252];
            XYZ = [Xn, Ylin, Zn];   % Nx3
            linRGB = (M * XYZ')';   % Nx3

            % Gamma correction (sRGB)
            rgb = color_lee.gamma_srgb(linRGB);

            % Clamp to [0,1]
            rgb = max(0, min(1, rgb));

            if n == 1
                rgb = rgb(:)'; % ensure 1x3 row vector
            end
        end

        function rgb = oklab(L, C, H)
            % oklab Converts OKLab LCH color to sRGB.
            % L: 0-1, C: chroma (~0-0.4), H: Hue 0-360 degrees
            % Returns Nx3 clamped sRGB values [0,1]

            if ~isscalar(L), L = L(:); end
            if ~isscalar(C), C = C(:); end
            if ~isscalar(H), H = H(:); end
            n = max([numel(L), numel(C), numel(H)]);
            if isscalar(L), L = repmat(L, n, 1); end
            if isscalar(C), C = repmat(C, n, 1); end
            if isscalar(H), H = repmat(H, n, 1); end

            % LCH -> OKLab (a,b)
            Hrad = H * pi / 180;
            a = C .* cos(Hrad);
            b = C .* sin(Hrad);

            % OKLab -> LMS (cube root space)
            M1 = [1,  0.3963377774,  0.2158037573;
                  1, -0.1055613458, -0.0638541728;
                  1, -0.0894841775, -1.2914855480];
            lab = [L, a, b]; % Nx3
            lms_ = (M1 * lab')'; % Nx3
            lms = lms_.^3;       % cube

            % LMS -> linear sRGB
            M2 = [ 4.0767416621, -3.3077115913,  0.2309699292;
                  -1.2684380046,  2.6097574011, -0.3413193965;
                  -0.0041960863, -0.7034186147,  1.7076147010];
            linRGB = (M2 * lms')'; % Nx3

            % Gamma + clamp
            rgb = color_lee.gamma_srgb(linRGB);
            rgb = max(0, min(1, rgb));

            if n == 1
                rgb = rgb(:)';
            end
        end

        function C_max = max_chroma(L, H)
            % max_chroma  Binary-search for the highest in-gamut chroma at given L, H.
            % L: 0-100, H: 0-360. Returns scalar or column vector.
            %
            % Example:
            %   C = color_lee.max_chroma(60, 140);
            %   rgb = color_lee.lch(60, C, 140);   % max-saturation green
            if ~isscalar(L), L = L(:); end
            if ~isscalar(H), H = H(:); end
            n = max(numel(L), numel(H));
            if isscalar(L), L = repmat(L, n, 1); end
            if isscalar(H), H = repmat(H, n, 1); end

            C_max = zeros(n, 1);
            for i = 1:n
                lo = 0; hi = 400;
                for iter = 1:40
                    mid = (lo + hi) / 2;
                    rgb = color_lee.lch_raw(L(i), mid, H(i));
                    if all(rgb >= -1e-6) && all(rgb <= 1 + 1e-6)
                        lo = mid;
                    else
                        hi = mid;
                    end
                end
                C_max(i) = lo;
            end
            if n == 1, C_max = C_max(1); end
        end

        function rgb = lch_maxchroma(L, H, fraction)
            % lch_maxchroma  Generate sRGB color at the gamut boundary.
            % L: 0-100, H: 0-360, fraction: 0-1 (default 1 = full boundary)
            %
            % Equal-lightness tricolor example:
            %   cred   = color_lee.lch_maxchroma(60,  10);
            %   cgreen = color_lee.lch_maxchroma(60, 140);
            %   cblue  = color_lee.lch_maxchroma(60, 260);
            if nargin < 3, fraction = 1.0; end
            C = color_lee.max_chroma(L, H) .* fraction;
            if isscalar(L) && ~isscalar(H)
                L = repmat(L, numel(H), 1);
            end
            rgb = color_lee.lch(L, C, H);
        end
    end

    methods (Static, Access = private)
        function out = gamma_srgb(linRGB)
            % sRGB gamma correction (linear -> gamma)
            mk = linRGB <= 0.0031308;
            out = 12.92 * linRGB .* mk + ...
                  (1.055 * max(linRGB, 0).^(1/2.4) - 0.055) .* (1 - mk);
        end

        function rgb = lch_raw(L, C, H)
            % lch_raw  Same as lch() but returns UNCLAMPED linear values.
            % Used internally by max_chroma binary search.
            Hrad = H * pi / 180;
            U = C * cos(Hrad);
            V = C * sin(Hrad);
            WP = [0.95047, 1.0, 1.08883];
            refd = WP(1) + 15*WP(2) + 3*WP(3);
            refU = 4*WP(1) / refd;
            refV = 9*WP(2) / refd;
            fY = (L + 16) / 116;
            ep = 216/24389; kp = 24389/27;
            Ylin = fY^3;
            if Ylin < ep, Ylin = (116*fY - 16) / kp; end
            if L == 0
                Un = refU; Vn = refV;
            else
                Un = U / (13*L) + refU;
                Vn = V / (13*L) + refV;
            end
            Xn = -(9 * Ylin * Un) / ((Un - 4) * Vn - Un * Vn);
            Zn = (9*Ylin - 15*Vn*Ylin - Vn*Xn) / (3*Vn);
            M = [3.2404542, -1.5371385, -0.4985314;
                -0.9692660,  1.8760108,  0.0415560;
                 0.0556434, -0.2040259,  1.0572252];
            rgb = (M * [Xn; Ylin; Zn])';
        end
    end
end
