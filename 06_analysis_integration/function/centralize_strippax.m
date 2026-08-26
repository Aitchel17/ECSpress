function out = centralize_strippax(pax_fwhm)
%CENTRALIZE_STRIPPAX  What centralized_paxfwhm keeps out of a line_fwhm.
%   IN   pax_fwhm  1x1 line_fwhm
%   OUT  out       1x1 struct

    % The one-dimensional series, and the two RAW kymographs.
    % Left behind: rotatecrop and mask, which are bulk nobody downstream reads,
    % and four of the six kymographs. Those four are regenerable from these two
    % plus the cutoffs already inside param, and one of them must not be used at
    % all -- kgph_*_processed clips each frame at its own 95th percentile, which
    % is narrower than the PVS crest, so the crest comes out flattened and the
    % renormalisation then calls that flat top 1.0. see CLAUDE_LOG.md
    keep = ["idx", "thickness", "displacement", "t_axis", "up_thicker", "param"];
    raw_kgph = ["kgph_pvs", "kgph_lumen"];

    out = struct();
    for f = keep
        if isprop(pax_fwhm, f) || isfield(pax_fwhm, f)
            out.(f) = pax_fwhm.(f);
        end
    end
    out.kymograph = struct();
    for f = raw_kgph
        if isfield(pax_fwhm.kymograph, f)
            out.kymograph.(f) = single(pax_fwhm.kymograph.(f));
        end
    end
end
