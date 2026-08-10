function [mask, info] = gate_rwidth(planes, u, v, min_snr)
%GATE_RWIDTH  Reject a vector SHORTER THAN THE PEAK that placed it. RETIRED.
%
% RETIRED 2026-08-06, taken out of analysis_pivensemble. It was already off by
% default long before that; this file only makes the comparison re-runnable.
%   err       a piece of tissue that did not move looks exactly like a short
%             vector, so the gate deletes the answer -- where the deformation
%             stops IS part of the result
%   see FINDINGS.md
%   why       in an ENSEMBLE the peak width carries the SPREAD of displacements
%             inside the window, not a localisation limit, so the ratio it forms
%             is not an SNR at all -- the name was wrong too
%
% IN   planes   corr_planes from piv_corr_ensemble; uses .maps
%      u, v     ny x nx float  the field as it stands
%      min_snr  float          reject where |d| < min_snr * peak half-width
% OUT  mask     ny x nx bool   true = reject
%      info     struct         on(bool) min_snr(float) peak_width(ny x nx float)

arguments
    planes  struct
    u       double
    v       double
    min_snr double = 1
end

info = struct('on', ~isempty(min_snr), 'min_snr', min_snr, 'peak_width', []);
mask = false(size(u));
if ~info.on; return; end

info.peak_width = piv_peak_width(planes.maps);
mask = hypot(u, v) < min_snr * info.peak_width;
end
