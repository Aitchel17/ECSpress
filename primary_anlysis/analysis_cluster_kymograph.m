function [cluster_column_idx,cluster_boundary]= analysis_cluster_kymograph(kymograph, numclusters)
% This function cla
% input: kymograph (vessel), number of cluster
% output: column order, cluster boundary

% 1. normalize preprocessing
zscore_kgph = (kymograph - mean(kymograph, 1)) ./ std(kymograph, [], 1);  % normalize each column
% 2. Distance claculation
D = pdist(zscore_kgph' , 'correlation'); % calculate pairwise correlation distance between all combinations of column
%% Clustering
tree_dist = linkage(D, 'average'); % Z: idx1, idx2, distance
cluster_idxs = cluster(tree_dist, 'maxclust', numclusters); % Cluster

%% Cluster sorting
% Using cluster mean values to sort clusters ( high mean --> dilated, low
% mean --> constricted)
cluster_means = zeros(numclusters, size(kymograph,1));
for k = 1:numclusters
    cluster_means(k,:) = mean(kymograph(:,cluster_idxs==k),2)';
end

% Inter cluster distance
cluster_dist = pdist(cluster_means, 'euclidean');
cluster_tree_dist = linkage(cluster_dist, 'average');
cluster_order = optimalleaforder(cluster_tree_dist, cluster_dist); 
% figure()
% dendrogram(cluster_tree_dist, reorder=cluster_order)
% Sort cluster
cluster_column_idx = [];
for k = 1:numel(cluster_order)
    idx_in_cluster = find(cluster_idxs==cluster_order(k));
    cluster_column_idx = [cluster_column_idx; idx_in_cluster];
end
% Cluster boundary
cluster_boundary = zeros(numclusters,2); 
current_idx = 1;

for k = 1:numel(cluster_order)
    idx_in_cluster = find(cluster_idxs==cluster_order(k));
    % disp(idx_in_cluster)
    num_elements = numel(idx_in_cluster);
    cluster_boundary(k,1) = current_idx;
    cluster_boundary(k,2) = current_idx + num_elements - 1;
    current_idx = current_idx + num_elements;
end

% %%
% kg_sorted = kymograph(:,cluster_column_idx);
% figure;
% subplot(1,2,1);
% imagesc(kymograph); colormap(turbo);
% title('Original Kymograph');
% axis tight;
% 
% subplot(1,2,2);
% imagesc(kg_sorted); colormap(turbo);
% title('PCA Sorted Kymograph');
% axis tight;
% % (위 코드 이어서 사용)
% 
% % kg_sorted와 final_order가 이미 만들어졌다고 가정하고 실행
% 
% % 클러스터별 경계 확인
% 
% 
% %% 시각화 (클러스터 경계선 추가)
% figure;
% imagesc(kg_sorted); colormap(turbo);
% title('Sorted Kymograph with Cluster Boundaries');
% axis tight; hold on;
% 
% % 클러스터 경계 표시 (세로선)
% for k = 1:size(cluster_boundary,1)-1
%     xline(cluster_boundary(k,2)+0.5,'m','LineWidth',1);
% end
% 
% % 클러스터 번호 추가
% for k = 1:size(cluster_boundary,1)
%     mid_point = mean(cluster_boundary(k,:));
%     text(mid_point, 5, sprintf('C%d', k), 'Color','white','FontSize',9,...
%         'HorizontalAlignment','center','FontWeight','bold');
% end
