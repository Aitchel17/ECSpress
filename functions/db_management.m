%%
dbfile = fullfile(pwd,"experiment_meta.db");
conn = sqlite(dbfile);
%%
stmt_ctable_mouse = strcat("CREATE TABLE", "IF NOT EXISTS", "Mouse","Mouse_ID TEXT PRIMARY KEY");
%%
backup_dir = uigetdir('','a');
%%
dir_info = dir(fullfile(backup_dir,'**','*.mdf'));
%%
dir_names = string({dir_info.folder}');
full_paths = fullfile(backup_dir,dir_names);
reg_pattern = "(?<Date>\d{6})_(?i:hql)(?<MouseNum>\d{2,3})(?:_?(?<Condition>.*))?";
%%
% 결과 변수 초기화
experiment_date = cell(size(dir_names));
mouse_id = cell(size(dir_names));
condition = cell(size(dir_names));
% 파싱 루프
for i = 1:numel(dir_names)
    tokens = regexp(dir_names{i}, reg_pattern, 'names');
    if ~isempty(tokens)
        experiment_date{i} = tokens.Date;
        mouse_id{i} = tokens.MouseNum;
        
        if isfield(tokens, 'Condition')
            condition{i} = tokens.Condition;
        else
            condition{i} = '';
        end
    else
        experiment_date{i} = '';
        mouse_id{i} = '';
        condition{i} = '';
    end
end

% 테이블로 정리
T = table(experiment_date, mouse_id, condition, dir_names, full_paths, ...
    'VariableNames', {'Date', 'MouseNum', 'Condition', 'FolderName', 'FullPath'});

%%
execute(conn, [
    "CREATE TABLE IF NOT EXISTS Mouse (" + ...
    "Mouse_ID TEXT PRIMARY KEY," + ...
    "Genotype TEXT," + ...
    "Notes TEXT" + ...
    ")"]);
%%