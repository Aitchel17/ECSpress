function thickness = fwhm_thickness(idx)
% Recover thickness (bv, uppvs, downpvs, totalpvs, eps) from cleaned upper and lower bv, pvs boundary

%   Caller: line_fwhm.getdiameter, tablegeneration_fwhmrelation 

%   IN   idx        1x1 struct   clean_upperBVboundary, clean_lowerBVboundary,
%                                clean_pvsupedge_idx, clean_pvsdownedge_idx
%   OUT  thickness  1x1 struct   bv, uppvs, downpvs, totalpvs, eps

    upper_bv = idx.clean_upperBVboundary;
    lower_bv = idx.clean_lowerBVboundary;
    up_pvs = idx.clean_pvsupedge_idx;
    down_pvs = idx.clean_pvsdownedge_idx;

    thickness.bv = lower_bv - upper_bv;
    thickness.uppvs = upper_bv - up_pvs;
    thickness.downpvs = down_pvs - lower_bv;
    thickness.totalpvs = thickness.uppvs + thickness.downpvs;
    thickness.eps = thickness.totalpvs + thickness.bv;
end
