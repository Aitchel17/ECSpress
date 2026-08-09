function [mask, info] = gate_rheight(planes, u, v, corr_thr)
%GATE_RHEIGHT  Reject a window whose correlation PEAK IS TOO LOW. RETIRED.
%
% RETIRED 2026-08-06, taken out of analysis_pivensemble. Kept because the
% comparison should be re-runnable, not re-argued.
%   see FINDINGS.md
%   why gone  gate_tomasi separates the same labels at no cost in LIVE windows,
%             and gate_corr_floor covers the messy-plane case without the bias
%
% IN   planes    corr_planes from piv_corr_ensemble; uses .corr
%      u, v      ny x nx float  the field as it stands
%      corr_thr  float          corr2 below this is rejected
% OUT  mask      ny x nx bool   true = reject
%      info      struct         on(bool) corr_thr(float) corr(ny x nx float)

arguments
    planes   struct
    u        double
    v        double
    corr_thr double = 0.5
end

info = struct('on', ~isempty(corr_thr), 'corr_thr', corr_thr, 'corr', []);
mask = false(size(u));
if ~info.on; return; end

info.corr = planes.corr;
% PIVlab's own filter, not a comparison of ours, so the threshold means what it
% means everywhere else in PIVlab
[filt_u, filt_v] = postproc.PIVlab_correlation_filter(u, v, corr_thr, info.corr);
mask = (isnan(filt_u) | isnan(filt_v)) & ~(isnan(u) | isnan(v));
end
