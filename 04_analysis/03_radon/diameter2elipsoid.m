function result = diameter2elipsoid(Dmat,theta)
%FIT2ELIPSOID Summary of this function goes here
%   Detailed explanation goes here
      

% INPUT
%   theta : [1 x Nangle] or [Nangle x 1] angles in radians, preferably 0~pi
%   Dmat  : [Nangle x Ntime] diameter data
%
% OUTPUT
%   result : struct with fields
%       .a           [Ntime x 1] semi-major axis
%       .b           [Ntime x 1] semi-minor axis
%       .phi         [Ntime x 1] orientation angle (rad)
%       .Dfit        [Ntime x Nangle] fitted diameter curves
%       .resnorm     [Ntime x 1] sum of squared residuals
%       .rmse        [Ntime x 1] root mean squared error
%       .exitflag    [Ntime x 1] optimizer exit flag
%       .AI_ratio    [Ntime x 1] b/a
%       .AI_1minus   [Ntime x 1] 1 - b/a
%       .AI_diffsum  [Ntime x 1] (a-b)/(a+b)

    [Ntheta, Nt] = size(Dmat);



    %% Output initialization
    a       = nan(Nt,1);
    b       = nan(Nt,1);
    phi     = nan(Nt,1);
    resnorm = nan(Nt,1);
    rmse    = nan(Nt,1);
    exitflag= nan(Nt,1);
    Dfit    = nan(Nt,Ntheta);

    % Optimization options
    opts = optimoptions('lsqnonlin', ...
        'Display', 'off', ...
        'MaxFunctionEvaluations', 5000, ...
        'MaxIterations', 500);

    for t = 1:Nt
        diameter1d = Dmat(:,t).';
        if numel(diameter1d) < 5
            continue;
        end

        % Initial guess from min/max
        dmax0 = max(diameter1d);
        dmin0 = min(diameter1d);

        a0 = dmax0 / 2;
        b0 = dmin0 / 2;

        % crude phi guess from angle at max diameter
        [~, idxMax] = max(diameter1d);
        phi0 = theta(idxMax);

        % enforce a >= b
        if b0 > a0
            tmp = a0; 
            a0 = b0; 
            b0 = tmp;
        end

        % ellipse orientation is pi-periodic
        phi0 = mod(phi0, pi);
        if phi0 >= pi/2
            phi0 = phi0 - pi;
        end
        p0 = [a0, b0, phi0];
        % bounds
        lb = [0,   0,   -pi/2];
        ub = [Inf, Inf,  pi/2];

        % robust scale based on MAD
        sigma = 1.4826 * mad(diameter1d, 1);
        if sigma <= 0 || ~isfinite(sigma)
            sigma = max(std(diameter1d), 1e-6);
        end
        robust_const = 2.5 * sigma;   % robust tuning constant

        % objective = robust-weighted residual
        objfun = @(p) robust_residual(p, theta, diameter1d, robust_const);
        % bounds
        lb = [eps, eps, -pi/2];
        ub = [Inf, Inf,  pi/2];
        try
            [pfit, ~, residual, ef] = lsqnonlin(objfun, p0, lb, ub, opts);

            % sort so a >= b
            aa = pfit(1);
            bb = pfit(2);
            pp = pfit(3);

            if bb > aa
                tmp = aa;
                aa = bb;
                bb = tmp;
                pp = pp + pi/2;
            end

            pp = mod(pp, pi);
            if pp >= pi/2
                pp = pp - pi;
            end
            yhat = ellipse_diameter_model([aa bb pp], theta);

            a(t) = aa;
            b(t) = bb;
            phi(t) = pp;
            Dfit(t,:) = yhat(:).';

            raw_res = diameter1d(:) - yhat(:);
            resnorm(t) = sum(raw_res.^2);
            rmse(t) = sqrt(mean(raw_res.^2));
            exitflag(t) = ef;

        catch
            % leave NaN if fit fails
        end
    end

    result = struct();
    result.a = a;
    result.b = b;
    result.phi = phi;
    result.Dfit = Dfit;
    result.resnorm = resnorm;
    result.rmse = rmse;
    result.exitflag = exitflag;

    result.AI_ratio   = b ./ a;           % circularity-like, closer to 1 = more circular
    result.AI_1minus  = 1 - b ./ a;       % 0 = circle
    result.AI_diffsum = (a - b) ./ (a + b); % 0 = circle
end


function r = robust_residual(p, theta, actual, c)
    a = p(1);
    b = p(2);
    phi = p(3);

    x = cos(theta - phi);
    y = sin(theta - phi);

    D = 2 * a * b ./ sqrt((b^2) .* x.^2 + (a^2) .* y.^2);
    error_term = actual - D;

    % pseudo-Huber-like weighting through residual transformation
    % lsqnonlin minimizes sum(r.^2), so use sqrt(weight).*e
    w = 1 ./ sqrt(1 + (error_term ./ c).^2);
    r = sqrt(w) .* error_term;
end





