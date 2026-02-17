function parsed_struct = mdf_commentparser(info_directory)
%MDF_COMMENTPARSER Extracts metadata from MDF comment strings
%   Parses comments for specific patterns:
%   - rcv/gcv (PMT gains): "rcv .53", "gcv.57", "rcv 50"
%   - Power: "20mW", "50mw"
%   - Depth: "30um", "89um"
%   - Zoom: "zoom 2.5", "x2.5" (if present)
%   - Vessel ID: First word if not a keyword, or explicit "PV01" etc.
% I would recommand to use AI tools to make regular expression, if
% additional information is in comment.

% Load MDF info

mdf = mdfExtractLoader(info_directory);

if isfield(mdf.info, 'Comments')
    comment_str = mdf.info.Comments;
else
    comment_str = '';
    warning('No Comments field found in MDF info.');
end

parsed_struct.comments = comment_str;

% Initialize fields
parsed_struct.vesselID = '';
parsed_struct.rPMT = [];
parsed_struct.gPMT = [];
parsed_struct.depth = [];
%
parsed_struct.power = [];
parsed_struct.power_percent = [];
parsed_struct.wavelength = [];

parsed_struct.zoom = [];
parsed_struct.objx = [];
parsed_struct.objy = [];
parsed_struct.objz = [];
parsed_struct.resolution = [];

if isempty(comment_str)
    return;
end

%% PMT Gain (Red/Green)
% Patterns: "rcv .53", "rcv.53", "rcv 0.53" -> 0.53V
% Case insensitive

% Red PMT (rcv)
r_match = regexp(comment_str, 'rcv\s*([.]?\d+(\.\d+)?)', 'tokens', 'ignorecase');
if ~isempty(r_match)
    val_str = r_match{1}{1};
    val = str2double(val_str);
    % Normalize: if > 1 (e.g. 5 -> 0.5, 47 -> 0.47)
    while val > 1
        val = val / 10;
    end
    parsed_struct.rPMT = strcat(string(val),' V');
end

% Green PMT (gcv)
g_match = regexp(comment_str, 'gcv\s*([.]?\d+(\.\d+)?)', 'tokens', 'ignorecase');
if ~isempty(g_match)
    val_str = g_match{1}{1};
    val = str2double(val_str);
    % Normalize: if > 1 (e.g. 5 -> 0.5, 47 -> 0.47)
    while val > 1
        val = val / 10;
    end
    parsed_struct.gPMT = strcat(string(val),' V');
end

%% Depth (um)
% Patterns: "30um", "89um", "50 um", "50.5um"
d_match = regexp(comment_str, '(\d+(\.\d+)?)\s*um', 'tokens', 'ignorecase');
if ~isempty(d_match)
    parsed_struct.depth = strcat(d_match{1}{1},' um');
end

%% Power (mW)
% Patterns: "20mW", "50mw", "20.5mW"
p_match = regexp(comment_str, '(\d+(\.\d+)?)\s*mW', 'tokens', 'ignorecase');
if ~isempty(p_match)
    parsed_struct.power = strcat(p_match{1}{1}, ' mW');
end
%% Vessel ID logic
% Detects: PV01 (Penetrating Vein), PA01 (Penetrating Arteriole),
%          SV01 (Surface Vein), SA01 (Surface Artery),
%          CA01 (Capillary), BP01 (Brain Parenchyma)
% Defaults to 'unknown' if not found.

parsed_struct.vesselID = 'unknown';

% Regex for allowed codes followed by digits
v_match = regexp(comment_str, '(PV|PA|SV|SA|CA|BP)\d+', 'match', 'ignorecase');

if ~isempty(v_match)
    % Take the first match and standardize to uppercase
    parsed_struct.vesselID = upper(v_match{1});
end
%% Zoom
parsed_struct.zoom = mdf.info.zoom;
%% resolution
parsed_struct.resolution = mdf.info.objpix;
%% Objective lens position
parsed_struct.objx = mdf.info.objx;
parsed_struct.objy = mdf.info.objy;
parsed_struct.objz = mdf.info.objz;

parsed_struct.power_percent = mdf.info.laserpower;
parsed_struct.wavelength = mdf.info.excitation;


end
