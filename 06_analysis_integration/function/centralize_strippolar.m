function out = centralize_strippolar(polarcluster)
%CENTRALIZE_STRIPPOLAR  What centralized_polarcluster keeps out of a polarcluster.
%   IN   polarcluster  1x1 struct
%   OUT  out           1x1 struct

    % Everything except clust_med_bv and clust_med_csf, which are the per-cluster
    % median image stacks and hold almost the whole file. cluster_bvcsf_constdil
    % is kept: it is the same kind of thing but only the constricted and dilated
    % pair, and the contour work reads it.
    keep = ["polar_profiles", "polar_theta", "cluster_summary", ...
        "cluster_num", "clusterboundary", "filtered_clusteridx", ...
        "constricted_cluster_idx", "dilated_cluster_idx", ...
        "constricted_medianimg", "dilated_medianimg", "cluster_bvcsf_constdil"];

    out = struct();
    for f = keep
        if isfield(polarcluster, f)
            out.(f) = polarcluster.(f);
        end
    end
end
