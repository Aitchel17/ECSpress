function result = heatmap_postprocessing(heatdata)
%HEATMAP_POSTPROCESSING Summary of this function goes here
%   Detailed explanation goes here

% 1. Extract main population of heatmap
tmp.bwmask = heatdata.xy_counts > 0;
tmp.cc = bwconncomp(tmp.bwmask);
tmp.stats = regionprops(tmp.cc);
tmp.area = [tmp.stats.Area];
[~,tmp.maxareaidx] = max(tmp.area);
tmp.cc.PixelIdxList{tmp.maxareaidx};
tmp.bwmask = zeros(size(tmp.bwmask));
tmp.bwmask(tmp.cc.PixelIdxList{tmp.maxareaidx}) = 1;

result.xy_maskedcounts = heatdata.xy_counts;
result.xy_maskedcounts(~tmp.bwmask) = NaN; % 마스크 적용
% Dev purpose, check data
% figure()
% imshow(heatdata.xy_maskedcounts)

% 2. NaN removal
result.xy_counts_clean = result.xy_maskedcounts;
% NaN removal for x axis
tmp.allnan_x = all(isnan(result.xy_counts_clean),1);
result.xy_counts_clean = result.xy_counts_clean(:,~tmp.allnan_x);
result.x_centers_clean = heatdata.x_centers(~tmp.allnan_x);
% NaN removal for y axis
tmp.allnan_y = all(isnan(result.xy_counts_clean),2);
result.xy_counts_clean = result.xy_counts_clean(~tmp.allnan_y,:);
result.y_centers_clean = heatdata.y_centers(~tmp.allnan_y);

% 3. Origin correction (Set baseline vessel, and PVS as origin)
% BV origin correction
tmp.baseline_bv = sum(result.xy_counts_clean,1,'omitmissing');
[~,tmp.maxbvloc]= max(tmp.baseline_bv);
result.x_baseceneters = result.x_centers_clean -  result.x_centers_clean(tmp.maxbvloc);
% pvs origin correction
tmp.baseline_pvs = sum(result.xy_counts_clean,2,'omitmissing');
[~,tmp.maxpvsloc]= max(tmp.baseline_pvs);
result.y_baseceneters = result.y_centers_clean - result.y_centers_clean(tmp.maxpvsloc);


% Mode of PVS for each x point
[~,result.modepvslocs] = max(result.xy_counts_clean,[],1);
result.modepvs = zeros([length(result.modepvslocs),1]);
    for idx = 1:length(result.modepvslocs)
        result.modepvs(idx) = result.y_baseceneters(result.modepvslocs(idx));
    end
end

