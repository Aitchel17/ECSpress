function io_saveinfo(info,save_folder)
    saveinfo = info;
    saveinfo.motionvertices = strjoin(string(reshape(saveinfo.motionvertices',1,[])));
    saveinfo.savedate = string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    
    infoFields = fieldnames(saveinfo);
    infoValues = struct2cell(saveinfo);
    table_info = table(infoFields, infoValues, 'VariableNames', {'Field', 'Value'});

    % Construct full file path
    save_infopath = fullfile(save_folder, [info.mdfName(1:end-4),'_info.txt']);
    % Write the table to an Excel file (overwrite the file initially)
    writetable(table_info, save_infopath);
end

