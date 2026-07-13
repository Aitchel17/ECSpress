function analysis_polar_makefig(polarcluster, save_dir)
% ANALYSIS_POLAR_MAKEFIG Renders the PAX-aligned polar contour profiles.
%   Plots the pre-computed profiles from analysis_clusterpolar_polarplot.
%   No transform here: it only iterates polarcluster.polar_profiles and plots.
%
%   Inputs:
%       polarcluster: Struct with .polar_theta (1x360) and .polar_profiles (1x4)
%       save_dir: Output directory (skip save if empty)

if ~isfield(polarcluster, 'polar_profiles') || isempty(polarcluster.polar_profiles)
    warning('polarcluster.polar_profiles missing. Run analysis_clusterpolar_polarplot first.');
    return;
end

theta = polarcluster.polar_theta;
% Rotate back to the ORIGINAL IMAGE angular frame (undo analyze_polar's PAX centering)
if isfield(polarcluster, 'polar_binstart_deg')
    theta = theta + deg2rad(polarcluster.polar_binstart_deg);
end

% Presentation only: colors by condition order (CBV, CPVS, DBV, DPVS)
clee = color_lee;
plot_colors = {clee.clist.red, clee.clist.darkgreen, clee.clist.magenta, clee.clist.green};

% Initialize polar figure
fig = make_fig('cluster_polar_contours', 'polar');
fig.bring_fig();
fig.update_figsize([6, 6]);
hold(fig.ax, 'on');

% Thin loop: plot each stored profile. NaN bins stay NaN so plot_polar
% (polarplot) breaks the line at missing angular sectors.
for i = 1:numel(polarcluster.polar_profiles)
    binned_r = polarcluster.polar_profiles(i).binned_r;
    color = plot_colors{min(i, numel(plot_colors))};
    fig.plot_polar(theta, binned_r, color, 'x');
end

% Mark the PAX axis direction on the outer edge: radial arrow + black circle.
if isfield(polarcluster, 'pax_angle')
    ax   = fig.ax;
    rmax = ax.RLim(2);
    th   = deg2rad(mod(-polarcluster.pax_angle, 360));   % PAX dir in image-aligned frame
    polarplot(ax, [th th], [0 rmax], 'k-', 'LineWidth', 1.2);              % radial line
    dth = deg2rad(6);  dr = 0.08 * rmax;                                    % arrowhead barbs
    polarplot(ax, [th th+dth], [rmax rmax-dr], 'k-', 'LineWidth', 1.2);
    polarplot(ax, [th th-dth], [rmax rmax-dr], 'k-', 'LineWidth', 1.2);
    polarplot(ax, th, rmax, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6); % black circle
    polarplot(ax, [th+pi th+pi], [0 rmax], 'k:', 'LineWidth', 0.8);         % opposite end (PAX is an axis)
    title(ax, sprintf('Cluster Contours (PAX angle = %.1f^{\\circ})', polarcluster.pax_angle));
else
    title(fig.ax, 'Cluster Contours (PAX angle)');
end

if ~isempty(save_dir)
    fig.save2svg(save_dir);
end

end
