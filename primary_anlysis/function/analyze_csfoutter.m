function [idx, kymographmask] = analyze_csfoutter(kymograph, bv_upidx, bv_downidx, threshold, offset)
% calculate_csf_boundaries - calculates upper and lower CSF boundaries based on given inputs.
%
% Inputs:
%   kymograph          - normalized CSF intensity profile (2D array)
%   bv_upidx  - upper boundary indices from blood vessel
%   bv_downidx - bottom boundary indices from blood vessel
%   threshold          - intensity threshold for boundary detection
%   offset        - initial thresholded binary mask
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
    edgebv_upkymograph = cropkymograph;
    edgebv_upkymograph(row_idx_grid > bv_upidx) = NaN;
    edgebv_downkymograph = cropkymograph;
    edgebv_downkymograph(row_idx_grid < bv_downidx) = NaN;
    
    cla(ax)
    imagesc(ax,edgebv_upkymograph);
    hold on
    plot(bv_upidx,'r')
  
    % 2.2 Shadow center detection (offset region) Outer boundary detection (upper)
    % Protect CSF region from wrong shadow detection below PVS signal
    %% 2.2.1 upmax position
    up_maxoffset = prctile(edgebv_upkymograph, 75, 1);
    up_maxidx = edgebv_upkymograph>=up_maxoffset;
    up_maxidx = up_maxidx .* row_idx_grid;
    up_maxidx(up_maxidx == 0) = NaN; % 1.6 convert 0 to NaN for median
    up_maxidx = round(median(up_maxidx,1,"omitmissing")); % Median quarter idx
    cla(ax)
    imagesc(ax,edgebv_upkymograph);
    plot(up_maxidx,'r')
    
    %% 2.2.1 downmax position
    down_maxoffset = prctile(edgebv_downkymograph, 75, 1);
    down_maxidx = edgebv_downkymograph>=down_maxoffset;
    down_maxidx = down_maxidx .* row_idx_grid;
    down_maxidx(down_maxidx == 0) = NaN; % 1.6 convert 0 to NaN for median
    down_maxidx = round(median(down_maxidx,1,"omitmissing")); % Median quarter idx
    cla(ax)
    imagesc(ax,edgebv_downkymograph);
    plot(down_maxidx,'r')

    %% 2.2.2 shadow center (if its just gradient cut off 10% of remaining region)
    edgebv_upmaxkymograph = edgebv_upkymograph;
    edgebv_upmaxkymograph(row_idx_grid > up_maxidx) = NaN;
    up_minoffset = prctile(edgebv_upmaxkymograph, 10, 1);
    up_minidx = edgebv_upmaxkymograph<= up_minoffset;
    up_minidx = up_minidx .* row_idx_grid;
    up_minidx(up_minidx == 0) = NaN; % 1.6 convert 0 to NaN for median
    up_minidx = round(median(up_minidx,1,"omitmissing")); % Median quarter idx
    cla(ax)
    imagesc(ax,edgebv_upmaxkymograph);
    hold on
    plot(up_minidx,'r')
    %%
    cla(ax)
    imagesc(ax,kymograph);
    hold on
    plot(up_minidx,'r')
    plot(up_maxidx,'r')
    
    %% 2.2.2 shadow center (if its just gradient cut off 10% of remaining region)
    edgebv_downmaxkymograph = edgebv_downkymograph;
    edgebv_downmaxkymograph(row_idx_grid < down_maxidx) = NaN;
    down_minoffset = prctile(edgebv_downmaxkymograph, 10, 1);
    down_minidx = edgebv_downmaxkymograph<= down_minoffset;
    down_minidx = down_minidx .* row_idx_grid;
    down_minidx(down_minidx == 0) = NaN; % 1.6 convert 0 to NaN for median
    down_minidx = round(median(down_minidx,1,"omitmissing")); % Median quarter idx
    cla(ax)
    imagesc(ax,edgebv_downmaxkymograph);
    hold on
    plot(down_minidx,'r')
    %%
    cla(ax)
    imagesc(ax,kymograph);
    hold on
    plot(up_minidx,'r')
    plot(up_maxidx,'r')

    plot(down_minidx,'r')
    plot(down_maxidx,'r')

    %% up PVS processing
    % from the scope (upmax position to local minimum position)
    up_mask = edgebv_upmaxkymograph;
    up_mask(row_idx_grid<up_minidx) = NaN;
    up_mask = up_mask -min(up_mask,[],1);
    up_mask = up_mask./max(up_mask,[],1);
    imagesc(ax, up_mask)

    %%
    % find half maximum points
    upthresholded = up_mask>=threshold;
    imagesc(ax, upthresholded)
    upthresholded_rawidx = upthresholded.*row_idx_grid;
    upthresholded_rawidx(upthresholded_rawidx==0) = Inf;
    upthresholded_rawidx = min(upthresholded_rawidx,[],1);
    %%
    upthresholded = imfill(upthresholded,'holes');
    nhood = [0 1 0; 0 0 0; 0 1 0];
    upthresholded = imerode(upthresholded,nhood);
    imagesc(ax,upthresholded)
    %%
    upperupboundary_idx = upthresholded.*row_idx_grid;
    upperupboundary_idx(upperupboundary_idx==0) = Inf;
    upperupboundary_idx = min(upperupboundary_idx,[],1);
    upperupboundary_idx(upperupboundary_idx == 0) = upthresholded_rawidx(upperupboundary_idx == 0);
    cla(ax)
    imagesc(ax,kymograph);
    hold on
    plot(upperupboundary_idx,'r')
  %%
    % Bottom kymograph processing
    % from the scope (upmax position to local minimum position)
    down_mask = edgebv_downmaxkymograph;
    imagesc(ax, down_mask)
    hold on
    plot(down_minidx,'r')
    %%
    down_mask(row_idx_grid>down_minidx) = NaN;
    down_mask = down_mask - min(down_mask,[],1);
    down_mask = down_mask./max(down_mask,[],1);
    imagesc(ax, down_mask)
    %%

    %% find half maximum points
    downthresholded = down_mask>=threshold;
    imagesc(ax, downthresholded)

    %%
    downrawthr_idx = downthresholded.*row_idx_grid;
    downrawthr_idx = max(downrawthr_idx,[],1);
    downthresholded = imfill(downthresholded,'holes');
    downthresholded = imerode(downthresholded,nhood);
    imagesc(ax,downthresholded)
    lowerlowboundary_idx = downthresholded.*row_idx_grid;
    lowerlowboundary_idx = max(lowerlowboundary_idx,[],1);
    lowerlowboundary_idx(lowerlowboundary_idx == 0) = downrawthr_idx(lowerlowboundary_idx == 0);
    %%
    cla(ax)
    imagesc(ax,kymograph);
    hold on
    plot(lowerlowboundary_idx,'r')
    plot(upperupboundary_idx,'r')
    %%
    
    %% output
    idx.pvs_upmax = up_maxidx;
    idx.pvs_downmax = down_maxidx;
    idx.pvs_upshadeloc = up_minidx;
    idx.pvs_downshadeloc = down_minidx;
    idx.pvs_upperboundary = upperupboundary_idx;
    idx.pvs_lowerboundary = lowerlowboundary_idx;
    %%
    uplocarray2d = repmat(idx.pvs_upperboundary, [sz(1), 1]); % 8.2
    upline = false(sz);
    upline(row_idx_grid == uplocarray2d) = 1;

    lowlocarray2d = repmat(idx.pvs_lowerboundary, [sz(1), 1]); % 8.2
    downline = false(sz);
    downline(row_idx_grid == lowlocarray2d) = 1;
    %%
    kymographmask.pvs_upline = upline;
    kymographmask.pvs_downline = downline;
    %%
    kymographmask.pvs_up = edgebv_upkymograph; % upper csf kymograph 
    kymographmask.pvs_up(row_idx_grid<upperupboundary_idx) = NaN;
    kymographmask.pvs_down = edgebv_downkymograph; % lower csf kymograph 
    kymographmask.pvs_down(row_idx_grid>lowerlowboundary_idx) = NaN; % lower csf kymograph 
    %%
end
