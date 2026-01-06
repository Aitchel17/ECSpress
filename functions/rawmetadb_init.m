function rawmetadb_init(backup_dir, dbfile)
% BUILD_MDF_DATABASE - 정리된 .mdf 실험 메타데이터를 sqlite DB로 저장하는 함수
% I believe when we recording with Sutter mScan, we usually make today's
% folder if not it will not work for you.
% The folder naming should be dddddd_cccddd_?????... corresponding to
% (YYMMDD_mouseid_experimenttype)
% 
% INPUTS:
%   backup_dir  - 최상위 실험 데이터 폴더 경로
%   dbfile      - 저장할 sqlite DB 파일 경로
%
% Requires: mdf_init
% 1. Check operating system in case db shifted to mainframe and access through my mac
% 2. Double check drive label and hdd label to be same (Check with your eye)
% 3. Gather directory of all .mdf file >>  
% 4. Make table
% 1
if ispc
    [drive_letter, ~, ~] = fileparts(backup_dir);
    [~, cmdout] = system(['vol ', drive_letter]);
    tokens = regexp(cmdout, 'Volume in drive [A-Z]: is (.*)\r?\n', 'tokens');
    if ~isempty(tokens)
        backup_name = strtrim(tokens{1}{1});
    else
        backup_name = 'UnknownDrive';
    end
elseif ismac || isunix
    parts = strsplit(backup_dir, '/');
    backup_name = parts{3};
else
    error('You are using strange OS');
end
% 2
if strcmpi(backup_name, 'UnknownDrive')
    backup_name = input('Drive name is not detected: ', 's');
else
    prompt = sprintf('Is this correct drive name? (I would not use default name) (y/n):', backup_name);
    confirm = input(prompt, 's');
    if lower(confirm) ~= 'y'
        backup_name = input('Type correct drive name (and change drive label after this): ', 's');
    end
end

% 3
all_mdf = dir(fullfile(backup_dir, '**', '*.mdf'));
folder_list = unique({all_mdf.folder});
Date = strings(size(folder_list));
MouseNum = strings(size(folder_list));
Condition = strings(size(folder_list));
PrimaryPath = folder_list';
PrimaryDrive = repmat(string(backup_name), size(folder_list));
reg_pattern = "(?<Date>\d{6})_(?i:hql)(?<MouseNum>\d{2,3})(?:_?(?<Condition>.*))?";
for i = 1:numel(folder_list)
    tokens = regexp(folder_list{i}, reg_pattern, 'names');
    if ~isempty(tokens)
        Date(i) = tokens.Date;
        MouseNum(i) = tokens.MouseNum;
        if isfield(tokens, 'Condition')
            Condition(i) = tokens.Condition;
        end
    end
end

FolderTable = table(Date, MouseNum, Condition, PrimaryPath, PrimaryDrive, 'VariableNames', {'Date','MouseNum','Condition','PrimaryPath','PrimaryDrive'});
%% XY Movie, Image Stack, Analog 테이블 분리 생성
XYMovieTable = table();
ImageStackTable = table();
AnalogTable = table();

file_names = {all_mdf.name}';
file_paths = fullfile({all_mdf.folder}', file_names);
[~, parent_idx] = ismember({all_mdf.folder}', folder_list);

for i = 1:numel(file_names)
    try
        [info, mobj] = mdf_init(all_mdf(i).folder, all_mdf(i).name);

        if strcmp(info.scanmode, 'XY Movie')
            beh_enable = strcmp(mobj.ReadParameter('Video Enabled'), '-1');

            xy_entry = table(i, string(file_names{i}), string(file_paths{i}), parent_idx(i), string(info.User), string(info.Date), beh_enable, string(info.fps), string(info.laserpower), string(info.fcount), string(info.fduration), string(backup_name), 'VariableNames', {'FileIdx','FileName','FilePath','ParentIdx','User','CreationDate','BehaviorEnable','FPS','LaserPower','FrameCount','FrameDuration','BackupDrive'});
            XYMovieTable = [XYMovieTable; xy_entry];

            % 아날로그 채널
            for ch = 0:7
                field_name = mobj.ReadParameter(sprintf('Analog Ch %d Name', ch));
                field_name(strfind(field_name,' ')) = '_';
                if ~isempty(field_name)
                    analog_entry = table(i, ch, string(field_name), string(mobj.ReadParameter(sprintf('Analog Ch %d Input Range', ch))), 'VariableNames', {'FileIdx', 'Channel', 'FieldName', 'InputRange'});
                    AnalogTable = [AnalogTable; analog_entry];
                end
            end

        elseif strcmp(info.scanmode, 'Image Stack')
            z_entry = table(i, string(file_names{i}), string(file_paths{i}), parent_idx(i), string(info.User), string(info.Date), string(mobj.ReadParameter('Average Count')),...
                               string(mobj.ReadParameter('Initial Intensity')),...
                               string(mobj.ReadParameter('Final Intensity')),...
                               string(mobj.ReadParameter('Z- interval')),...
                               string(backup_name), 'VariableNames', {'FileIdx','FileName','FilePath','ParentIdx','User','CreationDate','AverageCount','InitialIntensity','FinalIntensity','Zinterval','BackupDrive'});
            ImageStackTable = [ImageStackTable; z_entry];
        end

        mobj.release;
    catch
        warning('mdf 열기 실패: %s', file_paths{i});
    end
end

%% sqlite 저장
conn = sqlite(dbfile,'create');
existing_folders = sqlread(conn, 'FolderTable');
for i = 1:height(FolderTable)
    match_idx = ismember(existing_folders(:, {'Date', 'MouseNum', 'Condition'}), FolderTable(i, {'Date', 'MouseNum', 'Condition'}));
    if ~any(match_idx)
        sqlwrite(conn, 'FolderTable', FolderTable(i,:));
        folder_id = sqlread(conn, 'SELECT last_insert_rowid()');
    else
        folder_id = existing_folders.FolderID(find(match_idx,1));
        BackupEntry = table(folder_id, FolderTable.PrimaryPath(i), FolderTable.PrimaryDrive(i), 'VariableNames', {'FolderID', 'BackupPath', 'BackupDrive'});
        sqlwrite(conn, 'BackupTable', BackupEntry);
    end
end

sqlwrite(conn, 'XYMovieTable', XYMovieTable);
sqlwrite(conn, 'ImageStackTable', ImageStackTable);
sqlwrite(conn, 'AnalogTable', AnalogTable);

close(conn);
fprintf('DB 저장 완료: %s\n', dbfile);
end

