function out = centralize_strip(product, content)
%CENTRALIZE_STRIP  What centralized_<product>.mat keeps out of the loaded file.
%   Caller: centralize_mat.read_product
%   IN   product  1x1 str      the .mat stem, as it appears in centralize_mat
%        content  1x1 object   what came out of the session's <product>.mat
%   OUT  out      1x1 struct | the content itself
    switch product
        case "paxfwhm"
            out = centralize_strippax(content);
        case "polarcluster"
            out = centralize_strippolar(content);
        case "analysis_analog"
            out = centralize_stripanalog(content);
        otherwise
            out = content;
    end
end

function out = centralize_strippax(pax_fwhm)
%CENTRALIZE_STRIPPAX  What centralized_paxfwhm keeps out of a line_fwhm.
%   IN   pax_fwhm  1x1 line_fwhm
%   OUT  out       1x1 struct

    % 1. name the series and the two raw kymographs; rotatecrop, mask and four kymographs stay behind
    % kgph_*_processed is never kept : its per-frame 95th-percentile clip flattens the PVS crest
    keep = ["idx", "thickness", "displacement", "t_axis", "up_thicker", "param"];
    raw_kgph = ["kgph_pvs", "kgph_lumen"];

    % 2. copy the series
    out = struct();
    for f = keep
        if isprop(pax_fwhm, f) || isfield(pax_fwhm, f)
            out.(f) = pax_fwhm.(f);
        end
    end
    % 3. copy the raw kymographs as single
    out.kymograph = struct();
    for f = raw_kgph
        if isfield(pax_fwhm.kymograph, f)
            out.kymograph.(f) = single(pax_fwhm.kymograph.(f));
        end
    end
end

function out = centralize_strippolar(polarcluster)
%CENTRALIZE_STRIPPOLAR  What centralized_polarcluster keeps out of a polarcluster.
%   IN   polarcluster  1x1 analysis_clusterpolar | 1x1 struct (a file written before the class)
%   OUT  out           1x1 struct

    % 1. turn a file written before the class into the class
    if isstruct(polarcluster)
        polarcluster = analysis_clusterpolar.from_struct(polarcluster);
    end
    % 2. name all but clust_med_bv and clust_med_csf, the per-cluster median stacks that are most of the file
    keep = ["polar_profiles", "polar_theta", "cluster_summary", ...
        "clusterboundary", "filtered_clusteridx", ...
        "constricted_cluster_idx", "dilated_cluster_idx", ...
        "intermediate_cluster_idx", ...
        "constricted_medianimg", "dilated_medianimg"];

    % 3. copy cluster_num to the top level, as the rows before the class hold it, then the named fields
    out = struct();
    out.cluster_num = polarcluster.param.cluster_num;
    for f = keep
        out.(f) = polarcluster.(f);
    end
end

function out = centralize_stripanalog(primary_analog)
%CENTRALIZE_STRIPANALOG  What centralized_analysis_analog keeps out of an analysis_analog.
%   A plain struct, so the centralized file loads without 02_othersignal's classdef.
%   IN   primary_analog  1x1 analysis_analog
%   OUT  out             1x1 struct   air (puff sessions only), ball, ecog, emg, ecogspectrum

    % 1. keep the air channel at its own rate, one bool per sample, only where a puff was delivered;
    %    centralize_bout reads its presence as a whisker session
    out = struct();
    if isfield(primary_analog.rawdata, 'raw_Air_puff1')
        air_on = primary_analog.rawdata.raw_Air_puff1 > 0;
        if any(air_on)
            out.air = struct('on', air_on, 'fs', primary_analog.samplingfrequency);
        end
    end
    % 2. keep the resampled ball, ecog and emg; the raw-rate traces stay behind
    out.ball = pick_fields(primary_analog, 'ball', ...
        ["ds_fps", "rs_taxis", "resampled_ball"]);
    out.ecog = pick_fields(primary_analog, 'ecog', ...
        ["ds_fps", "resampled_ecog", "resampled_taxis"]);
    out.emg = pick_fields(primary_analog, 'emg', ...
        ["ds_fps", "rs_taxis", "resampled_power", "resampled_signal"]);

    % 3. keep the spectrogram axes and one power array
    out.ecogspectrum = struct();
    if isprop(primary_analog, 'ecog') && isfield(primary_analog.ecog, 'ecogspectrum')
        spectrum = primary_analog.ecog.ecogspectrum;
        for f = ["T", "F", "params", "movingwin"]
            if isfield(spectrum, f)
                out.ecogspectrum.(f) = spectrum.(f);
            end
        end
        % baseline_S (sleep: log ratio to awake, F x T) or else S (whisker: raw power, T x F), each under its own name
        if isfield(spectrum, 'baseline_S')
            out.ecogspectrum.baseline_S = single(spectrum.baseline_S);
        elseif isfield(spectrum, 'S')
            out.ecogspectrum.S = single(spectrum.S);
        end
    end
end

function out = pick_fields(container, group, keep)
%PICK_FIELDS  The fields of container.(group) named in keep that exist; empty when the group is missing.
    out = struct();
    if ~isprop(container, group) || isempty(container.(group))
        return
    end
    for f = keep
        if isfield(container.(group), f)
            out.(f) = container.(group).(f);
        end
    end
end
