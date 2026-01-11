function polar_boundary = polar_halfmax(sinogram_complex,bvthr,bvoffset,pvsthr,pvsoffset)
%POLAR_HALFMAX Summary of this function goes here
%   Detailed explanation goes here
polar_boundary = struct();
const_bvsinogram = squeeze(sinogram_complex(:,:,1,1));
dilate_bvsinogram = squeeze(sinogram_complex(:,:,1,2));
% follow similar strategy of fwhm based analysis of bottom side
polar_boundary.const_bvboundary = getbvboundary(const_bvsinogram,bvthr,bvoffset);
polar_boundary.dilate_bvboundary = getbvboundary(dilate_bvsinogram,bvthr,bvoffset);
% PVS boundary calculation
const_pvsinogram = squeeze(sinogram_complex(:,:,2,1));
dilate_pvsinogram = squeeze(sinogram_complex(:,:,2,2));
polar_boundary.const_pvsboundary = getpvsboundary(const_pvsinogram, polar_boundary.const_bvboundary,pvsthr,pvsoffset);
polar_boundary.dilate_pvsboundary = getpvsboundary(dilate_pvsinogram, polar_boundary.dilate_bvboundary,pvsthr,pvsoffset);

%% PVS analysis using igkl channel
% Debug purpose
% figure()
% %
% imagesc(sinogram)
% colormap gray
% %
% hold on
% plot(boundary_idx, 'r')
% %
% plot(downoffsetloc,'y')
% %
% plot(bvboundary,'y')

end

function boundary_idx = getbvboundary(sinogram, bvthr, bvoffset)
    sz = size(sinogram);
    [row_idx_grid, ~] = ndgrid(1:sz(1), 1:sz(2));
    down_offset = prctile(sinogram,bvoffset,1);
    downoffsetloc =  sinogram<=down_offset;
    downoffsetloc = row_idx_grid.*downoffsetloc;
    downoffsetloc(downoffsetloc == 0) = Inf;
    downoffsetloc = min(downoffsetloc,[],1);
    downoffsetloc = medfilt1(downoffsetloc,3);
    sinogram(row_idx_grid > downoffsetloc) = NaN;
    sinogram = sinogram - down_offset;
    sinogram = sinogram./max(sinogram,[],1);
    sinogram_thr = sinogram > bvthr;
    boundary_idx = row_idx_grid.*sinogram_thr; % 7.2
    boundary_idx = max(boundary_idx, [], 1); % 7.4
    boundary_idx = squeeze(boundary_idx);
end

function boundary_idx = getpvsboundary(sinogram,vesselboundary,pvsthr,pvsoffset)
    sz = size(sinogram);
    [row_idx_grid, ~] = ndgrid(1:sz(1), 1:sz(2));
    sinogram(row_idx_grid<vesselboundary) = NaN;
        
    down_offset = prctile(sinogram,pvsoffset,1);
    downoffsetloc =  sinogram<=down_offset;
    downoffsetloc = row_idx_grid.*downoffsetloc;
    downoffsetloc(downoffsetloc == 0) = Inf;
    downoffsetloc = min(downoffsetloc,[],1);
    downoffsetloc = medfilt1(downoffsetloc,3);
    %
    sinogram(row_idx_grid > downoffsetloc) = NaN;
    sinogram = sinogram - down_offset;
    sinogram = sinogram./max(sinogram,[],1);
    %
    sinogram_thr = sinogram > pvsthr;
    boundary_idx = row_idx_grid.*sinogram_thr; % 7.2
    boundary_idx = max(boundary_idx, [], 1); % 7.4
    boundary_idx = squeeze(boundary_idx);
    % If there is no contrast between PVS and ISF, drop the angle
        % downoffsetloc: 
    drop_criteria = downoffsetloc - vesselboundary;
    drop_criteria = drop_criteria < 3;
    boundary_idx(drop_criteria) = NaN;
end

