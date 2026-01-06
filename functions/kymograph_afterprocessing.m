function output = kymograph_afterprocessing(linedata_session)
%KYMOGRAPH_AFTERPROCESSING Summary of this function goes here
%   Detailed explanation goes here
disp(linedata_session.infodict('Comments'))
%
input.pax = linedata_session.fwhmline.pax;
input.kgph_bv = input.pax.kymograph.kgph_bv;
input.kgph_csf = input.pax.kymograph.kgph_csf;
%1 using ceterbvidx, seperate bottom half and upperhalf
%2 using minimum position identification might screwed <<< lets try
%without min position 
input.centerbvidx = ceil((input.pax.idx.bv_lowerboundary+input.pax.idx.bv_upperboundary)/2);
input.rowidx = input.pax.mask.rowidx;

% median filter
tmp.upbv = medfilt1(input.pax.idx.bv_upperboundary,20);
tmp.upcell = medfilt1(input.pax.idx.pvs_upshadeloc,20);
tmp.uppvs = medfilt1(input.pax.idx.pvs_upperboundary,20);
tmp.downbv = medfilt1(input.pax.idx.bv_lowerboundary,20);
tmp.downcell = medfilt1(input.pax.idx.pvs_downshadeloc,20);
tmp.downpvs = medfilt1(input.pax.idx.pvs_lowerboundary,20);

% remove first point to avoid medfilt
tmp.downbv = tmp.downbv(2:end);
tmp.downcell = tmp.downcell(2:end);
tmp.downpvs = tmp.downpvs(2:end);
tmp.upbv = tmp.upbv(2:end);
tmp.upcell = tmp.upcell(2:end);
tmp.uppvs = tmp.uppvs(2:end);

% thicknesscalculation
output.bv_thickness = tmp.downbv - tmp.upbv;
output.bv_changes = output.bv_thickness-min(output.bv_thickness);

output.uppvs_thickness = tmp.upbv - tmp.uppvs;
output.uppvs_changes = output.uppvs_thickness-min(output.uppvs_thickness);

output.downpvs_thickness = tmp.downpvs - tmp.downbv;
output.downpvs_changes = output.downpvs_thickness-min(output.downpvs_thickness);

output.totalpvs_changes = output.uppvs_changes + output.downpvs_changes;
output.totalpvs_thickness = output.downpvs_thickness+output.uppvs_thickness;
end

