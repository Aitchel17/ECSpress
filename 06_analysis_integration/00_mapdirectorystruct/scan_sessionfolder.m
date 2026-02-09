function session_data = scan_sessionfolder(session_dir)
% SCAN_SESSIONFOLDER Scan a session folder for files and metadata
%
% Input:
%   session_dir - Full path to session directory
%
% Output:
%   session_data - Struct containing parsed metadata, main files, and analysis sub-structs

session_data = struct();
session_data.directory = session_dir;

if ~isfolder(session_dir)
    warning('Session folder not found: %s', session_dir);
    session_data.error = true;
    return;
end

%% 1. Parse Metadata (MDF parser)
parsed_meta = mdf_commentparser(session_dir);
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
