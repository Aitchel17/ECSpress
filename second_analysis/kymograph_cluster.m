%% analysis window initialization
window1=figure(Name='window1',NumberTitle='off');
window2=figure(Name='window2',NumberTitle='off');
window3=figure(Name='window3',NumberTitle='off');
window4=figure(Name='window4',NumberTitle='off');

%% gaussian fit full width half max
test.snippet_fname = 'HQL072_250331_003';
test.pax = linedata_struct.(test.snippet_fname).fwhmline.pax;
test.kgph_bv = test.pax.kymograph.kgph_bv;
test.kgph_csf = test.pax.kymograph.kgph_csf;
%1 using ceterbvidx, seperate bottom half and upperhalf
%2 using minimum position identification might screwed <<< lets try
%without min position 
test.centerbvidx = ceil((test.pax.idx.bv_lowerboundary+test.pax.idx.bv_upperboundary)/2);
test.rowidx = test.pax.mask.rowidx;
%%
test.kgph_csf_up = test.kgph_csf;
test.kgph_csf_up(test.rowidx<=test.centerbvidx) = NaN;
test.kgph_csf_down = test.kgph_csf;
test.kgph_csf_down(test.rowidx>=test.centerbvidx) = NaN;
%%
figure(window1)
imagesc(test.kgph_csf_up)

figure(window2)
imagesc(test.kgph_csf_down)
%%
figure(window3)
imagesc(test.kgph_bv)
%%

tmp.norm_kgph_csf = test.kgph_csf./max(test.kgph_csf,[],1);
%%
% 입력: kymograph 이미지 (예: 70 x 7500)
kg = test.kgph_bv;  % size: height x time

% 1. 전처리: 정규화 (각 column마다)
kg_noedxrm = (kg - mean(kg, 1)) ./ std(kg, [], 1);  % normalize each column

% 2. 거리 계산 (column 간 유사도)
dist = pdist(kg_norm', 'euclidean');  % transpose: compare columns
Z = linkage(dist, 'average');

%% 3. 최적 순서 계산 (비슷한 column끼리 정렬)
order = optimalleaforder(Z, dist);

% 4. 정렬된 kymograph 재구성
kg_sorted = kg(:, order);  % 시간 축 정렬

% 5. 시각화
figure;
subplot(1,2,1);
imagesc(kg); title('Original'); colormap(turbo); axis tight;

subplot(1,2,2);
imagesc(kg_sorted); title('Sorted by column similarity'); colormap(turbo); axis tight;
%%
kg =  test.kgph_csf;  % size: height x time

kg_norm = (kg - mean(kg, 1)) ./ std(kg, [], 1);  % normalize each column

[~, score] = pca(gather(kg_gpu)); % GPU에서 PCA는 MATLAB에서 제한적이라 gather필요

[~, order] = sort(score(:,1)); % 가장 강력한 패턴 기준 정렬

kg_sorted = kg(:,order);

figure;
subplot(1,2,1);
imagesc(kg); colormap(turbo);
title('Original Kymograph');
axis tight;

subplot(1,2,2);
imagesc(kg_sorted); colormap(turbo);
title('PCA Sorted Kymograph');
axis tight;
%%
%% 입력 kymograph
kg = test.kgph_bv; % height x time (e.g., 70 x 7500)

% GPU 기반 정규화 및 거리계산
kg_gpu = gpuArray(single(kg')); % (time x height)

% 정규화: 각 column을 평균 0, 표준편차 1로 조정
kg_gpu = (kg_gpu - mean(kg_gpu,2))./std(kg_gpu,[],2);

% pdist2를 이용한 빠른 거리 행렬 계산
D_gpu = pdist2(kg_gpu, kg_gpu, 'euclidean');
D = gather(D_gpu);
%%
kg = test.kgph_bv; % height x time (e.g., 70 x 7500)
% column vector smoothing
kg_smooth = movmean(kg, 5, 1); % window=5, 각 column 방향 smooth
% 정규화
kg_norm = (kg_smooth - mean(kg_smooth, 1)) ./ std(kg_smooth, [], 1);

% correlation distance: check correlation to match form
D = pdist(kg_norm' , 'correlation');
% Clustering
maxClusters = 2; % num of cluster
Z = linkage(D, 'average');
T = cluster(Z, 'maxclust', maxClusters);
% 클러스터 대표값 계산 및 클러스터 간 재정렬
cluster_means = zeros(maxClusters, size(kg,1));
for k = 1:maxClusters
    cluster_means(k,:) = mean(kg(:,T==k),2)';
end

% 클러스터 간 거리 행렬 계산 후 정렬
cluster_dist = pdist(cluster_means, 'euclidean');
Z_cluster = linkage(cluster_dist, 'average');

% dendrogram 대신 leaves 순서 직접 얻기
cluster_order = optimalleaforder(Z_cluster, cluster_dist);

% 최종 column 정렬 순서 만들기 (클러스터 수에 맞춰 반복)
final_order = [];
numActualClusters = numel(cluster_order); % 실제 클러스터 수에 맞춤
for k = 1:numActualClusters
    idx_in_cluster = find(T==cluster_order(k));
    final_order = [final_order; idx_in_cluster];
end

%
kg_sorted = kg(:,final_order);

figure;
subplot(1,2,1);
imagesc(kg); colormap(turbo);
title('Original Kymograph');
axis tight;

subplot(1,2,2);
imagesc(kg_sorted); colormap(turbo);
title('PCA Sorted Kymograph');
axis tight;
% (위 코드 이어서 사용)

% kg_sorted와 final_order가 이미 만들어졌다고 가정하고 실행

% 클러스터별 경계 확인
cluster_boundary = zeros(maxClusters,2); % 각 클러스터의 시작과 끝 인덱스
current_idx = 1;

for k = 1:numel(cluster_order)
    idx_in_cluster = find(T==cluster_order(k));
    num_elements = numel(idx_in_cluster);

    cluster_boundary(k,1) = current_idx;
    cluster_boundary(k,2) = current_idx + num_elements - 1;

    current_idx = current_idx + num_elements;
end

%% 시각화 (클러스터 경계선 추가)
figure;
imagesc(kg_sorted); colormap(turbo);
title('Sorted Kymograph with Cluster Boundaries');
axis tight; hold on;

% 클러스터 경계 표시 (세로선)
for k = 1:size(cluster_boundary,1)-1
    xline(cluster_boundary(k,2)+0.5,'m','LineWidth',1);
end

% 클러스터 번호 추가
for k = 1:size(cluster_boundary,1)
    mid_point = mean(cluster_boundary(k,:));
    text(mid_point, 5, sprintf('C%d', k), 'Color','white','FontSize',9,...
        'HorizontalAlignment','center','FontWeight','bold');
end



%%
figure(window4)
plot(kg_sorted(:,1000:1500))
%%
figure(window4)
plot(mean(kg_sorted(:,800:1000),2))


%%
tmp.centerbvidx(1)

%%
test_pax.idx.pvs_downshadeloc()

%%
pax = test_pax.mask.pvs_up;
tmp.downcsfkgph = test_pax.mask.pvs_down;
%%
tmp.downcsfkgph = test_pax.kymograph.kgph_csf;
%%


%%

figure()

for idx = 1: 1000

plot(tmp.downcsfkgph(tmp.centerbvidx:test_pax.idx.pvs_downshadeloc,idx))
pause(0.1)

end

%%
for ixd =1:214
plot(tmp.upcsfkgph(:,idx))
hold on
plot(tmp.downcsfkgph(:,idx))
pause(0.5)
hold off
end


