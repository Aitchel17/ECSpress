function [idx, kymographmask] = get_pvsoutter(kymograph, bv_upidx, bv_downidx)
% calculate_csf_boundaries - calculates upper and lower CSF boundaries based on given inputs.
%
% Inputs:
%   kymograph          - normalized CSF intensity profile (2D array)
%   bv_upidx  - upper boundary indices from blood vessel
%   bv_downidx - bottom boundary indices from blood vessel
%
% Outputs:
%   idx.up_csf_boundary    - upper CSF boundary indices
%   idx.down_csf_boundary  - lower CSF boundary indices
%   kymographmask.cb_thresholded     - updated binary mask with CSF boundaries

%%
fig = figure('Name','Check kymograph'); % 0.1.1 initialize figure
ax = axes(fig); % 0.1.2 initialize axis
imagesc(ax,kymograph) % 0.1.3 show original kymograph for interactive crop
%%
% 0. output struct
idx = []; % 1D row index per slice
kymographmask = []; % 2D Mask corresponding to kymograph or masked kymograph

% 1. size and grid matrix generation
sz = size(kymograph);
[row_idx_grid, ~] = ndgrid(1:sz(1), 1:sz(2)); % 2D row index grid
% bv_centeridx = ceil((bv_downidx+bv_upidx)/2); % center point to
% divide kymograph  251013. CSF boundary restraint

%% Restrict scope

while true
    cla(ax)
    imagesc(ax,kymograph)
    upcropval = input('Get upper crop idx: '); % 0.2.1 Get input
    upcropkymograph = kymograph;

    upcropkymograph(1:upcropval,:) = NaN;
    imagesc(ax,upcropkymograph)
    flag = input('looks ok? type "y":  ', 's'); % 0.3.5 check
    if strcmp(flag,'y') % 0.3.5 check
        break
    end
end

while true
    cla(ax)
    imagesc(ax,kymograph)
    downcropval = input('Get lower crop idx: '); % 0.2.1 Get input
    cropkymograph = upcropkymograph;
    cropkymograph(downcropval:end,:) = NaN;
    imagesc(ax,cropkymograph)
    flag = input('looks ok? type "y":  ', 's'); % 0.3.5 check
    if strcmp(flag,'y') % 0.3.5 check
        break
    end
end



%% 2.Separation process (CSF external boundary using vessel upper boundary)
% 2.1 Upper vessel kymograph generation (above upper vessel boundary)
bv_mididx = round((bv_upidx+bv_downidx)/2);

% upper vessel edge side kymograph
edgebv_upkymograph = cropkymograph;
edgebv_upkymograph(row_idx_grid > bv_upidx) = NaN; % upper vessel edge kymograph
up_maxidx = findmax(edgebv_upkymograph); %
edgebv_upmaxkymograph = edgebv_upkymograph;
edgebv_upmaxkymograph(row_idx_grid > up_maxidx) = NaN;
[upedge_confined,upedge_minidx] = findshadow(edgebv_upmaxkymograph,true); % 3. minidx
[upedge_thresholded, upedge_rawidx] = halfmax_confinedkymograph(upedge_confined,true); % upcenter, boundary
[upedge_idx,~] = clean_thresholdedkymograph(upedge_thresholded>0, true, ax);
%%
figure(fig)
imagesc(ax,upedge_confined)
%%


%%
% upper vessel center side kymograph
centbv_upkymograph = cropkymograph; % upper vessel center kymograph
centbv_upkymograph(row_idx_grid > bv_mididx) = NaN;
centbv_upmaxkymograph = centbv_upkymograph;
centbv_upmaxkymograph(row_idx_grid < up_maxidx) = NaN;
[upcent_thresholded, upcent_rawidx] = halfmax_confinedkymograph(centbv_upmaxkymograph,false); % upcenter, boundary
[upcent_idx,clean_kymograph] = clean_thresholdedkymograph(upcent_thresholded>0, false, ax);
%% debuging purpose (center bv vs edge bv scope setup)
figure(fig)
imagesc(ax,upedge_thresholded>0)
hold on
plot(bv_upidx,'r')
plot(bv_mididx,'g')
%%
imagesc(ax,upedge_confined)
hold on
plot(up_maxidx, 'r')
plot(upedge_minidx, 'g')
plot(bv_mididx, 'g')
%% Check scope of up-center
imagesc(ax,centbv_upmaxkymograph)
hold on
plot(up_maxidx, 'r')
plot(upedge_minidx, 'g')
plot(bv_mididx, 'g')
%%
imagesc(ax,cropkymograph)
hold on
plot(upedge_rawidx, 'r')
plot(upcent_rawidx, 'g')
%%
imagesc(ax,cropkymograph)
hold on
plot(upedge_idx, 'r')
plot(upcent_idx, 'g')
%% lower vessel edge side kymograph
edgebv_downkymograph = cropkymograph; % lower vessel edge kymograph
edgebv_downkymograph(row_idx_grid < bv_downidx) = NaN;
down_maxidx = findmax(edgebv_downkymograph);
edgebv_downmaxkymograph = edgebv_downkymograph; % crop to maximum
edgebv_downmaxkymograph(row_idx_grid < down_maxidx) = NaN;
[downedge_confined,downedge_minidx] = findshadow(edgebv_downmaxkymograph,false);
[downedge_thresholded, downedge_rawidx] = halfmax_confinedkymograph(downedge_confined,false); % downcenter, boundary
%%

sum(downedge_thresholded,1);


%%
figure(fig)

imagesc(ax,cropkymograph)

hold on
plot(ax,bv_downidx,'g')
%%
plot(ax,down_maxidx,'r')

%%
[downedge_idx,~] = clean_thresholdedkymograph(downedge_thresholded, false, ax);

% lower vessel center side kymograph
centbv_downkymograph = cropkymograph;
centbv_downkymograph(row_idx_grid < bv_mididx) = NaN;
centbv_downmaxkymograph = centbv_downkymograph;
centbv_downmaxkymograph(row_idx_grid > down_maxidx) = NaN;
[downcent_thresholded, downcent_rawidx] = halfmax_confinedkymograph(centbv_downmaxkymograph,true); % downcenter, boundary
[downcent_idx,clean_kymograph] = clean_thresholdedkymograph(downcent_thresholded, true, ax);
%%
%% debuging purpose (center bv vs edge bv scope setup)
figure(fig)
imagesc(ax,cropkymograph)
hold on
plot(bv_downidx,'r')
plot(bv_mididx,'g')
%%
cla(ax)
imagesc(ax,downedge_confined)
hold on
plot(down_maxidx, 'r')
plot(downedge_minidx, 'g')
plot(bv_mididx, 'g')
%% Check scope of down-center
cla(ax)
imagesc(ax,centbv_downmaxkymograph)
hold on
plot(down_maxidx, 'r')
plot(downedge_minidx, 'g')
plot(bv_mididx, 'g')

%%
imagesc(ax,cropkymograph)
hold on
plot(downedge_rawidx, 'r')
plot(downcent_rawidx, 'g')
%%
imagesc(ax,cropkymograph)
hold on
plot(downedge_idx, 'r')
plot(downcent_idx, 'g')
%%
imagesc(ax,cropkymograph)
hold on
plot(upedge_idx, 'r')
plot(upcent_idx, 'g')
plot(downedge_idx, 'r')
plot(downcent_idx, 'g')
%%
imagesc(ax,cropkymograph)
hold on
plot(upedge_idx, 'r')
plot(upcent_idx, 'g')
plot(downedge_idx, 'r')
plot(downcent_idx, 'g')



%% output
idx.pvsupedge_idx = upedge_idx;
idx.pvsupcent_idx = upcent_idx;
idx.pvsdownedge_idx = downedge_idx;
idx.pvsdowncent_idx = downcent_idx;

idx.pvsupedge_rawidx   = upedge_rawidx;
idx.pvsupcent_rawidx   = upcent_rawidx;
idx.pvsdownedge_rawidx = downedge_rawidx;
idx.pvsdowncent_rawidx = downcent_rawidx;

%%
uplocarray2d = repmat(idx.pvsupedge_idx, [sz(1), 1]); % 8.2
upline = false(sz);
upline(row_idx_grid == uplocarray2d) = 1;

lowlocarray2d = repmat(idx.pvsdownedge_idx, [sz(1), 1]); % 8.2
downline = false(sz);
downline(row_idx_grid == lowlocarray2d) = 1;
%%
kymographmask.pvs_upline = upline;
kymographmask.pvs_downline = downline;
%%
imagesc(ax,kymographmask.pvs_upline+kymographmask.pvs_downline)
%%
kymographmask.pvs_up = edgebv_upkymograph; % upper csf kymograph
kymographmask.pvs_up(row_idx_grid<upedge_idx) = NaN;
kymographmask.pvs_down = edgebv_downkymograph; % lower csf kymograph
kymographmask.pvs_down(row_idx_grid>downedge_rawidx) = NaN; % lower csf kymograph

imagesc(ax,sum(cat(3,kymographmask.pvs_up,kymographmask.pvs_down),3,'omitnan'));

end

function [boundary_idx,clean_kymograph] = clean_thresholdedkymograph(thresholded, upmod, debug_ax)
x_length = size(thresholded,2);
clean_kymograph = imfill(thresholded,'holes');
boundary_idx = zeros([1,x_length]);

nhood = zeros([7,3]);
if upmod
    nhood(1:2,2) =1; % up
else
    nhood(end-1:end,2) =1; % down
end
clean_kymograph = imerode(clean_kymograph,nhood);
clean_kymograph = imdilate(clean_kymograph,nhood);
clean_kymograph(1:2,:) = 0; % remove up dilation
clean_kymograph(end-1:end,:) = 0; % remove bottom dilation

%%
zero_column = sum(clean_kymograph,1);
zero_column = zero_column ==0;
%%
clean_kymograph(:,zero_column) = thresholded(:,zero_column);
%%
zeros_array = sum(clean_kymograph,1);
%%

%%
%%
cla(debug_ax);
%plot(thresholded(:,8))
% %%
imagesc(debug_ax,clean_kymograph)
%%
for c_idx = 1: x_length
    szfilt.x = clean_kymograph(:,c_idx);
    szfilt.dx =diff([0;szfilt.x;0]);
    szfilt.sx = find(szfilt.dx == 1);
    szfilt.ex = find(szfilt.dx == -1);
    szfilt.lx = szfilt.ex -szfilt.sx;

    szfilt.y = ~szfilt.x;
    szfilt.dy =diff([0;szfilt.y;0]);
    szfilt.sy = find(szfilt.dy == 1);
    szfilt.ey = find(szfilt.dy == -1);
    szfilt.ly = szfilt.ey -szfilt.sy;
    [~,szfilt.lxmaxszidx] = max(szfilt.lx);
    if upmod
        szfilt.count = 1;% up
        szfilt.ly = szfilt.ly(2:end); % up
        while szfilt.count < szfilt.lxmaxszidx
            if szfilt.lx(szfilt.count) < szfilt.ly(szfilt.count)
                % disp(szfilt.lx(szfilt.count))
                szfilt.sx(1:szfilt.count) = Inf;
                szfilt.ex(1:szfilt.count) = Inf;
            end
            szfilt.count=szfilt.count+1;
        end
        boundary_idx(c_idx) = min(szfilt.sx);
    else
        szfilt.count = numel(szfilt.lx); % down
        szfilt.ly = szfilt.ly(1:end-1); % down
        while szfilt.count > szfilt.lxmaxszidx
            if szfilt.lx(szfilt.count) < szfilt.ly(szfilt.count)
                %disp(szfilt.lx(szfilt.count))
                szfilt.sx(szfilt.count:end) = 0;
                szfilt.ex(szfilt.count:end) = 0;
            end
            szfilt.count=szfilt.count-1;
        end
        boundary_idx(c_idx) = max(szfilt.ex);
    end

end
%%
%%
cla(debug_ax);
%%
imagesc(debug_ax,thresholded)
%%
imagesc(debug_ax,clean_kymograph)
%%
end

function maxidx = findmax(cropkymograph)
sz = size(cropkymograph);
[row_idx_grid, ~] = ndgrid(1:sz(1), 1:sz(2)); % 2D row index grid
maxoffset = prctile(cropkymograph, 75, 1);
maxidx = cropkymograph>=maxoffset;
maxidx = maxidx .* row_idx_grid;
maxidx(maxidx == 0) = NaN; % 1.6 convert 0 to NaN for median
maxidx = round(median(maxidx,1,"omitmissing")); % Median quarter idx
maxidx = medfilt1(maxidx,11);
end

function [crop_mask,minidx] = findshadow(maxcropkymograph,upmod)
% upmod true : minposition is above maximum, so remove upper portion become
% NaN
sz = size(maxcropkymograph);
[row_idx_grid, ~] = ndgrid(1:sz(1), 1:sz(2)); % 2D row index grid
minoffset = prctile(maxcropkymograph, 10, 1);
minidx = maxcropkymograph<= minoffset;
minidx = minidx .* row_idx_grid;
minidx(minidx == 0) = NaN; % 1.6 convert 0 to NaN for median
minidx = round(median(minidx,1,"omitmissing")); % Median quarter idx
minidx = medfilt1(minidx,11);
crop_mask = maxcropkymograph;
if upmod
    %%
    mincolumnidx = maxcropkymograph>0;
    mincolumnidx = mincolumnidx.*row_idx_grid;
    %%
    mincolumnidx(mincolumnidx ==0) = Inf;
    %%
    mincolumnidx = min(mincolumnidx,[],1);
    %%
    differential_column = mincolumnidx - minidx;
    overshadow_column = differential_column <5;
    minidx(overshadow_column) = mincolumnidx(overshadow_column); % if thickness thinner than 5 pixel, don't narrow down
    minidx(overshadow_column) = mincolumnidx(overshadow_column);

    %%
    crop_mask(row_idx_grid<minidx) = NaN;
else
    maxcolumnidx = maxcropkymograph>0;
    maxcolumnidx = maxcolumnidx.*row_idx_grid;
    maxcolumnidx = max(maxcolumnidx,[],1);
    differential_column = minidx - maxcolumnidx;
    overshadow_column = differential_column <5;
    minidx(overshadow_column) = maxcolumnidx(overshadow_column); % if thickness thinner than 5 pixel, don't narrow down
    minidx(overshadow_column) = maxcolumnidx(overshadow_column);
    crop_mask(row_idx_grid > minidx) = NaN;
end
%%

end

function [thresholded_kgph,thresholded_rawidx] = halfmax_confinedkymograph(confined_kymograph,upmod)
sz = size(confined_kymograph);
[row_idx_grid, ~] = ndgrid(1:sz(1), 1:sz(2)); % 2D row index grid
normconfined_kgph = confined_kymograph -min(confined_kymograph,[],1);
zerocolumn = sum(normconfined_kgph,1,'omitmissing');
zerocolumn = zerocolumn ==0 ;
normconfined_kgph(:,zerocolumn) = confined_kymograph(:,zerocolumn);
normconfined_kgph = normconfined_kgph./max(normconfined_kgph,[],1);
%%
zerocolumn = sum(normconfined_kgph,1,'omitmissing');
zerocolumn = zerocolumn ==0 ;
%%

thresholded_kgph = normconfined_kgph>=0.5;


thresholded_rawidx = thresholded_kgph.*row_idx_grid;
if upmod
    thresholded_rawidx(thresholded_rawidx==0) = Inf;
    thresholded_rawidx = min(thresholded_rawidx,[],1);
else
    thresholded_rawidx = max(thresholded_rawidx,[],1);
end
%%

end

function [sz1d, sz1dind,sizematrix] = binary2columnsize(binary_matrix)
sz = size(binary_matrix);
transit = abs(diff(binary_matrix,1,1));
transit = [true(1,sz(2));transit];
%%
columnszind = cumsum(transit,1);
%%
sz1dind = columnszind + (0:sz(2)-1)*max(columnszind,[],'all'); % make each column to be unique 1d count and reconstruction
sz1dind = sz1dind(:);
sz1d = accumarray(sz1dind, 1); % count idx = size
%%
sizematrix = sz1d(sz1dind); % duplicate size back to correspoinding idx
sizematrix = reshape(sizematrix,sz(1), sz(2));
end

