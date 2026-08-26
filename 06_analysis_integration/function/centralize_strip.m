function out = centralize_strip(product, content)
%CENTRALIZE_STRIP  What centralized_<product>.mat keeps out of the loaded file.
%   IN   product  1x1 str      the .mat stem, as it appears in centralize_primary
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
