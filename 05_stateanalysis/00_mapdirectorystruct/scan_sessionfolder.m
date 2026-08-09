function [session_data, parsed_meta] = scan_sessionfolder(session_dir, cached_meta)
% SCAN_SESSIONFOLDER Scan a session folder for files and metadata
%
% Input:
%   session_dir - Full path to session directory
%   cached_meta - (optional) previously parsed MDF metadata struct for this
%                 session. MDF acquisition metadata is immutable, so when a
%                 cached struct is supplied the expensive MDF parse is skipped.
%
% Output:
%   session_data - Struct containing parsed metadata, main files, and analysis sub-structs
%   parsed_meta  - The MDF metadata struct used (freshly parsed or the cached
%                  input), so the caller can persist it for reuse. Empty if the
%                  folder was missing.

if nargin < 2
    cached_meta = [];
end

session_data = struct();
session_data.directory = session_dir;
parsed_meta = [];

if ~isfolder(session_dir)
    warning('Session folder not found: %s', session_dir);
    session_data.error = true;
    return;
end

%% 1. Parse Metadata (MDF parser) -- immutable per session; reuse cache if given
if isempty(cached_meta)
    parsed_meta = mdf_commentparser(session_dir);
else
    parsed_meta = cached_meta;
end
fields = fieldnames(parsed_meta);
for i = 1:numel(fields)
    % Sanitize name -> if parser returns 'VesselID', use it
    session_data.(fields{i}) = string(parsed_meta.(fields{i}));
end


%% 2. Check Dictionary Contents (Main Files)
contents = dir(session_dir);
% Filter out . and ..
contents = contents(~ismember({contents.name}, {'.', '..'}));

% Expected suffixes/patterns for main files
% Map: FieldName -> Suffix/Pattern
main_files_map = {
    'Analog', '_analog.txt';
    'Ch1', '_ch1.tif';
    'Ch2', '_ch2.tif';
    'Eye', '_eye.avi';
    'Info', '_info.txt';
    'Motion', '_motion.txt';
    'Whisker', '_whisker.avi';
    'Peripheral', 'peripheral';
    'PrimaryAnalysis', 'primary_analysis'
    'State_analysis', 'state_analysis'
    };

for i = 1:size(main_files_map, 1)
    field = main_files_map{i, 1};
    pattern = main_files_map{i, 2};

    matches = contents(endsWith({contents.name}, pattern));
    if ~isempty(matches)
        session_data.(field) = matches(1).name;
    else
        session_data.(field) = 'NA';
    end
end

end
