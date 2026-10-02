function out = centralize_strippolar(polarcluster)
%CENTRALIZE_STRIPPOLAR  What centralized_polarcluster keeps out of a polarcluster.
%   IN   polarcluster  1x1 analysis_clusterpolar | 1x1 struct (a file written before the class)
%   OUT  out           1x1 struct

    % Everything except clust_med_bv and clust_med_csf, which are the per-cluster
    % median image stacks and hold almost the whole file. The two medianimg are
    % the same kind of thing for the constricted and dilated picks only.
    if isstruct(polarcluster)
        polarcluster = analysis_clusterpolar.from_struct(polarcluster);
    end
    keep = ["polar_profiles", "polar_theta", "cluster_summary", ...
        "clusterboundary", "filtered_clusteridx", ...
        "constricted_cluster_idx", "dilated_cluster_idx", ...
        "intermediate_cluster_idx", ...
        "constricted_medianimg", "dilated_medianimg"];

    % cluster_num at the top level, as in the rows written before the class
    out = struct();
    out.cluster_num = polarcluster.param.cluster_num;
    for f = keep
        out.(f) = polarcluster.(f);
    end
end
