function io_saveanalog(analog, save_folder, info)
    % Auto detection of saving function
    
    % Ball saving - Only if 'ds_ball' exists
    if isfield(analog, 'ds_ball')
        table_ball = array2table(analog.ds_ball', 'VariableNames', {'times (s)', 'Ball velocity (V)'}); 
        save_ballpath = fullfile(save_folder, [info.mdfName(1:end-4), '_ball.txt']);
        writetable(table_ball, save_ballpath);
    else
        disp('Warning: Field "ds_ball" is missing, skipping Ball data saving.');
    end
    
    % EMG saving - Only if 'ds_EMG' exists
    if isfield(analog, 'ds_EMG')
        table_emg = array2table(analog.ds_EMG', 'VariableNames', {'times (s)', 'EMG (V)'}); 
        save_emgpath = fullfile(save_folder, [info.mdfName(1:end-4), '_emg.txt']);
        writetable(table_emg, save_emgpath);
    else
        disp('Warning: Field "ds_EMG" is missing, skipping EMG data saving.');
    end
    
    % ECoG saving - Only if 'ecog_spectrum' exists
    if isfield(analog, 'ecog_spectrum')
        % Prepare ECoG data for saving
        time_spectrum = array2table(analog.ecog_spectrum.log_norm_spectrum);
        time_spectrum.Properties.VariableNames = string(analog.ecog_spectrum.t_axis); 
        f_column = array2table(flip(analog.ecog_spectrum.f_axis'));
        f_column.Properties.VariableNames{1} = 'Frequency(Hz)\Times (sec)';
        table_ecog = [f_column, time_spectrum];
        save_ecogpath = fullfile(save_folder, [info.mdfName(1:end-4), '_ecog.xlsx']);
        writetable(table_ecog, save_ecogpath, 'Sheet', 'ecog_mtspec', 'WriteVariableNames', true);
    else
        disp('Warning: Field "ecog_spectrum" is missing, skipping ECoG data saving.');
    end
    
    % Plotting and saving figure
    % Check if any of the relevant fields for plotting exist
    if isfield(analog, 'ds_ball') || isfield(analog, 'ds_EMG') || isfield(analog, 'ecog_spectrum')
        [analog_fig, ~] = analog_plot(analog);  % Generate the plot and get the figure handle

        % Save the figure in .fig format
        save_analogsummarypath = fullfile(save_folder, [info.mdfName(1:end-4), '_analog.fig']);
        savefig(analog_fig, save_analogsummarypath);

        % Save the figure as a vectorized .pdf
        save_analogsummarypath_pdf = fullfile(save_folder, [info.mdfName(1:end-4), '_analog.pdf']);
        exportgraphics(analog_fig, save_analogsummarypath_pdf);
    else
        disp('Warning: No valid fields found for plotting (ds_ball, ds_EMG, or ecog_spectrum), skipping plot saving.');
    end
end
