function [idx, kymograph_mask,parameter] = get_bvoutter(kymograph)
%% 0. Manual crop kymograph, initialize interactive window
% 0.0 default paramter
threshold = 0.5;
upoffset = 5;
lowoffset = 5;
% 0.1 initialize user interface
fig = figure('Name','Check kymograph'); % 0.1.1 initialize figure
ax = axes(fig); % 0.1.2 initialize axis
imagesc(ax,kymograph) % 0.1.3 show original kymograph for interactive crop
% 0.2 up crop idx determination
while true
    figure(fig)
    cla(ax) % 0.3.1 reset axis to prevent overlay burden
    imagesc(ax,kymograph) % 0.3.2 plot entire kymograph
    upcropval = input('Get upper crop idx: '); % 0.2.1 Get input
    imagesc(ax,kymograph(upcropval:end,:)) % 0.2.2 show up crop result
    upcropflag = input('looks ok? type "y":  ', 's'); % 0.2.3 Check
    if strcmp(upcropflag,'y') % 0.2.3
        break %0.1.5
    end
end
% 0.3 down crop idx determination
cla(ax) % 0.3.1 reset axis to prevent overlay burden
imagesc(ax,kymograph) % 0.3.2 plot entire kymograph
while true
    cla(ax) % 0.3.1 reset axis to prevent overlay burden
    imagesc(ax,kymograph) % 0.3.2 plot entire kymograph
    lowcropval = input('Get lower crop idx: '); % 0.3.3 get input
    imagesc(ax,kymograph(1:lowcropval,:)) % 0.3.4 show down crop result
    upcropflag = input('looks ok? type "y":  ', 's'); % 0.3.5 check
    if strcmp(upcropflag,'y') % 0.3.5 check
        break
    end
end
%% 0. result: cropped kymograph

cropkymograph = kymograph(upcropval:lowcropval,:); % 0.4
sz = size(cropkymograph); % 1.1 initialize parameters % 0.4
[row_idx_grid, ~] = ndgrid(1:sz(1), 1:sz(2)); % 0.4
max_value = max(cropkymograph,[],1);
%% 1. Estimate center of vessel
cla(ax) % 1.1 reset kymograph to show cropped kymograph
imagesc(ax,cropkymograph) % 1.2 show cropped kymograph
%%
% Max val calculation
max_offset = prctile(cropkymograph,75,1); % 1.3 get quarter max value
max_idx =  cropkymograph>=max_offset; % 1.4 threshold quarter max
max_idx = max_idx.*row_idx_grid; % 1.5 get idx
max_idx(max_idx == 0) = NaN; % 1.6 convert 0 to NaN for median
max_idx = round(median(max_idx,1,"omitmissing")); % Median quarter idx
filtered_max_idx = medfilt1(max_idx,31); % smoothing max position
hold on
plot(ax,filtered_max_idx,'r')


%% 2. Upkymograph analysis
norm_kymograph = cropkymograph./max_value;
upkymograph = norm_kymograph;
upkymograph(row_idx_grid > max_idx) = NaN;
downkymograph = norm_kymograph;
downkymograph(row_idx_grid < max_idx) = NaN;
%%
cla(ax) % 1.1 reset kymograph to show cropped kymograph
imagesc(ax,upkymograph) % 1.2 show cropped kymograph

while true
    up_offset = prctile(upkymograph,upoffset,1); %% offset
    upoffsetloc =  upkymograph<=up_offset;
    upoffsetloc = row_idx_grid.*upoffsetloc;
    upoffsetloc(upoffsetloc == 0) = NaN; % 1.6 convert 0 to NaN for median
    upoffsetloc = round(median(upoffsetloc,1,"omitmissing")); % Median quarter idx
    upoffsetloc = medfilt1(upoffsetloc,31);
    cla(ax) % 1.1 reset kymograph to show cropped kymograph
    imagesc(ax,upkymograph) % 1.2 show cropped kymograph
    hold on
    plot(upoffsetloc, 'r')
    upoffsetflag = input('looks ok? "y" ','s');
    if strcmp(upoffsetflag,'y')
        break
    else
        upoffset = input('Default offset was 5% put new offset (%)');
    end
end
% 7. bottom proessing


while true
    cla(ax) % 1.1 reset kymograph to show cropped kymograph
    imagesc(ax,downkymograph) % 1.2 show cropped kymograph
    down_offset = prctile(downkymograph,lowoffset,1);
    downoffsetloc =  downkymograph<=down_offset;
    downoffsetloc = row_idx_grid.*downoffsetloc;
    downoffsetloc(downoffsetloc == 0) = NaN; % 1.6 convert 0 to NaN for median
    downoffsetloc = round(median(downoffsetloc,1,"omitmissing")); % Median quarter idx
    downoffsetloc = medfilt1(downoffsetloc,31);
    hold on
    plot(downoffsetloc,'r-')
    flag = input('looks ok? "y" ', 's');
    if strcmp(flag,'y')
        break
    else
        lowoffset = input('Default offset was 5% put new offset (%)');
    end
end
%%
upkymograph(row_idx_grid < upoffsetloc) = NaN;
downkymograph(row_idx_grid > downoffsetloc) = NaN;

%%
cla(ax) % 1.1 reset kymograph to show cropped kymograph
imagesc(ax,upkymograph) % 1.2 show cropped kymograph
while true
    upkymograph = upkymograph - min(upkymograph,[],1); % below offset become negative
    upkymograph = upkymograph./max(upkymograph,[],1); % Max to be 1
    upkymograph_thr = upkymograph > threshold; % 6.2 thresholding
    upperboundary_idx = row_idx_grid .* (upkymograph_thr); % & row_idx_grid <= maxlocarray2d); % 6.3
    upperboundary_idx(upperboundary_idx==0) = Inf;
    upperboundary_idx = min(upperboundary_idx, [], 1); % 6.4
    upperboundary_idx = squeeze(upperboundary_idx);
    cla(ax) % 1.1 reset kymograph to show cropped kymograph
    imagesc(ax,upkymograph) % 1.2 show cropped kymograph
    hold on
    plot(upperboundary_idx, 'r')
    upboundaryflag = input('looks ok? "y" ', 's');
    if strcmp(upboundaryflag,'y')
        disp('continue to bottom processing')
    else
        threshold = input('Default threshold is 0.5 (FWHM) put new threshold');
    end

    cla(ax) % 1.1 reset kymograph to show cropped kymograph
    imagesc(ax,downkymograph) % 1.2 show cropped kymograph
    downkymograph = downkymograph - min(downkymograph,[],1);
    downkymograph = downkymograph./max(downkymograph,[],1);
    cla(ax) % 1.1 reset kymograph to show cropped kymograph
    imagesc(ax,downkymograph) % 1.2 show cropped kymograph
    downkymograph_thr = downkymograph>threshold; % 6.2
    sprintf('thresholding at %.2f',threshold)
    lowerboundary_idx = row_idx_grid.*downkymograph_thr; % 7.2
    lowerboundary_idx = max(lowerboundary_idx, [], 1); % 7.4
    lowerboundary_idx = squeeze(lowerboundary_idx);
    hold on
    plot(lowerboundary_idx,'r')
    flag = input('looks ok? "y" ', 's');
    if strcmp(flag,'y')
        break
    else
        threshold = input('Default threshold is 0.5 (FWHM) put new threshold');
    end
end


%% 8. Final mask processing

%%


% 9. return
parameter.upcropidx = upcropval;
parameter.lowcropidx = lowcropval;
parameter.threshold = threshold;
parameter.upoffset = upoffset;
parameter.lowoffset = lowoffset;

idx =struct();
idx.max_idx = filtered_max_idx + upcropval;
idx.upoffsetloc = upoffsetloc + upcropval;
idx.downoffsetloc = downoffsetloc + upcropval;
idx.upperBVboundary = upperboundary_idx + upcropval;
idx.lowerBVboundary = lowerboundary_idx + upcropval;
cla(ax,'reset')
imagesc(ax,kymograph)
hold on
plot(idx.upperBVboundary,'r')
plot(idx.lowerBVboundary, 'r')

ori_size = size(kymograph); % 1.1 initialize parameters % 0.4
[ori_row_idx_grid, ~] = ndgrid(1:sz(1), 1:sz(2)); % 0.4

ori_upkymograph = NaN(ori_size);
ori_upkymograph(upcropval:lowcropval,:) = upkymograph;
ori_downkymograph = NaN(ori_size);
ori_downkymograph(upcropval:lowcropval,:) = downkymograph;


uplocarray2d = repmat(idx.upperBVboundary, [sz(1), 1]); % 8.2
upline = false(sz);
upline(ori_row_idx_grid == uplocarray2d) = 1;

downlocarray2d = repmat(idx.lowerBVboundary, [sz(1), 1]); % 8.4
downline = false(sz);
downline(ori_row_idx_grid == downlocarray2d) = 1;

midline = false(sz);
midline(row_idx_grid == ceil((downlocarray2d+uplocarray2d)/2)) = 1;

kymograph_mask = struct();
kymograph_mask.rowidx = row_idx_grid;
kymograph_mask.up = ori_upkymograph;
kymograph_mask.down = ori_downkymograph;
kymograph_mask.upline = upline;
kymograph_mask.downline = downline;
kymograph_mask.midline = midline;
end

