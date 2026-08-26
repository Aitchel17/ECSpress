function thickness = fwhm_thickness(idx)
%FWHM_THICKNESS  The lumen and perivascular series, from the four boundary rows.
%   Every series below is a difference of two of the four rows line_fwhm cleans,
%   which is why centralized_paxfwhm_state carries the rows and not these.
%   line_fwhm.getdiameter calls this, so the arithmetic exists once and a reader
%   working from the centralized rows gets the same five series the class does.
%   IN   idx        1x1 struct   clean_upperBVboundary, clean_lowerBVboundary,
%                                clean_pvsupedge_idx, clean_pvsdownedge_idx
%   OUT  thickness  1x1 struct   bv, uppvs, downpvs, totalpvs, eps
%                                SAME shape and class as the rows it was given --
%                                this only subtracts, so line_fwhm gets back exactly
%                                what its five inline lines used to produce
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
