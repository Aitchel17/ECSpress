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
%% Acquisition metadata, one field at a time
% err      only zoom used to be guarded, so a session whose folder has no
%          *_info.txt died on objpix and took the whole mapdirstruct with it.
%          NONE of the folders that hit it has a paxfwhm.mat -- the scan was
%          dying on folders it had nothing to collect from
%          see FINDINGS.md
% note     the fields keep the [] they were initialised with when absent, so a
%          session with no metadata reads as empty rather than as a wrong number
info_fields = { ...
    'zoom',          'zoom'; ...
    'resolution',    'objpix'; ...
    'objx',          'objx'; ...
    'objy',          'objy'; ...
    'objz',          'objz'; ...
    'power_percent', 'laserpower'; ...
    'wavelength',    'excitation'};
for fi = 1:size(info_fields, 1)
    if isfield(mdf.info, info_fields{fi, 2})
        parsed_struct.(info_fields{fi, 1}) = mdf.info.(info_fields{fi, 2});
    end
end

if isempty(comment_str)
    return;
else

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
end


end
