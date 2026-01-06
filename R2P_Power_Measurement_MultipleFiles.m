%%  Resonant 2Photon Power (R2P) Measurement
%   Problem:
%       Measuring the power of the resonant 2photon by manually altering
%       the power on the MScan software and recording the signal from the
%       power meter is relatively inefficient (time intensive and subject
%       to examiner error). 
%       By setting Mscan to average multiple frames and using the automatic 
%       intensity control (linear) to alter the power intensity, a graph of
%       the R2P power intensity vs laser output can be created.
%   Directions:
%       Create a folder that contains all .txt files obtained from the
%       ThorLabs Power Meter. This script will prompt you to select that
%       folder.
%   Output:
%       Two main outputs: a set of graphs and a summary table. Each .txt
%       file will be plotted and specified input laser intensity
%       percentages will be shown with respect to the output laser
%       intensity. The summary table will provide data that can be copied
%       to the excel sheet on google drive to record the R2P.
%
function R2P_Power_Measurement_MultipleFiles


%% System Properties
SystemProperties.FilePath.GetFolder = uigetdir;                                                 % Choose folder path location
    SystemProperties.FilePath.Address = SystemProperties.FilePath.GetFolder + "\*.txt";         % Convert to filepath
    SystemProperties.FilePath.Folder = dir(SystemProperties.FilePath.Address);                  % Identify the folder directory 
    SystemProperties.FilePath.Data.Length = length(SystemProperties.FilePath.Folder);           % Determine the number of files in folder directory
    SystemProperties.FilePath.Data.Address = erase(SystemProperties.FilePath.Address,"*.txt");  % Create beginning address for file path
SystemProperties.Threshold = 1;                                                                 % Adjust minimum threshold to remove data before/after active recording

%% Data Processing
for i = 1:SystemProperties.FilePath.Data.Length                                                                                 % Begin for loop to iterate over individual files
    SystemProperties.FilePath.Data.Name{i} = erase(SystemProperties.FilePath.Folder(i).name,".txt");                            % Create names of individual files based only on MScan wavelength
    R2Pm.RawData{i} = readtable(fullfile(SystemProperties.FilePath.Data.Address, SystemProperties.FilePath.Folder(i).name));    % Read data from individual files (full file path for each file contained)
    
    R2Pm.Data{i} = table;                                                                               % Create table for processing data 
        R2Pm.Data{i}.DateTime = R2Pm.RawData{i}.Var1;                                                   % Insert time (DateTime format) variable into table
        R2Pm.Data{i}.ReferenceTime = datenum(R2Pm.RawData{i}.Var1);                                     % Insert reference time (time since recording starts) variable into table. Convert from datetime to number format
        R2Pm.Data{i}.Intensity = R2Pm.RawData{i}.Var2*1000;                                             % Insert laser intensity variable into table. Multiply x1000 to convert power to mW units
        R2Pm.Data{i} = R2Pm.Data{i}(~(R2Pm.Data{i}.Intensity < SystemProperties.Threshold),:);          % Remove all rows that give power intensity below threshold
            R2Pm.Data{i}.DateTime = R2Pm.Data{i}.DateTime - R2Pm.Data{i}.DateTime(1);                   % Reset time (DateTime format) variable to 0
            R2Pm.Data{i}.ReferenceTime = R2Pm.Data{i}.ReferenceTime - R2Pm.Data{i}.ReferenceTime(1);    % Reset time (number format) variable to 0
    R2Pm.Results{i} = table;                                                                                        % Create table for results
        [R2Pm.Results{i}.OutputIntensity, index] = findpeaks([R2Pm.Data{i}.Intensity]);                             % Find local maximum laser power intensities and associated index
        R2Pm.Results{i}.Time = R2Pm.Data{i}.ReferenceTime(index);                                                   % Determine time variable based on index
        R2Pm.Results{i}.InputIntensity = [0, 100*((1:size(R2Pm.Results{i},1)-1) / (size(R2Pm.Results{i},1)-1))]';   % Create input power intensity based on percentages (actually input for MScan)
    R2Pm.Analysis{i} = table;                                                                                                                                           % Create table for data analysis
        R2Pm.Analysis{i}.InterestValues = [0:5:70,80:10:100]';                                                                                                          % Power values of interest as defined by excel sheet
            for j = 1:length(R2Pm.Analysis{i}.InterestValues)                                                                                                           % Begin for loop to iterate over values of interest
                R2Pm.Analysis{i}.Time(j,1) = interp1(R2Pm.Results{i}.InputIntensity, R2Pm.Results{i}.Time,R2Pm.Analysis{i}.InterestValues(j));                          % Interpolate input intensity for scanning time based on values of interest
                R2Pm.Analysis{i}.OutputIntensity(j,1) = interp1(R2Pm.Results{i}.InputIntensity, R2Pm.Results{i}.OutputIntensity,R2Pm.Analysis{i}.InterestValues(j));    % Interpolate input intensity for output intensity based on values of interest 
            end                                                                                                                                                         % End for loop
    R2Pm.Final.Final = table;                                       % Create table of useful information
        R2Pm.Final.Raw(:,1) = R2Pm.Analysis{1}.InterestValues;      % First column provides the power input percentages
        R2Pm.Final.Raw(:,i+1) = R2Pm.Analysis{i}.OutputIntensity;   % All following columns provide the measured output intensity
        R2Pm.Final.Final = array2table(R2Pm.Final.Raw);             % Collect useful data into table
end                                                                                                                             % End for loop

%% Plot & Output

for i = 1:SystemProperties.FilePath.Data.Length
    figure(i)                                                                                                               % Create plots
    plot(R2Pm.Data{i}.DateTime, R2Pm.Data{i}.Intensity, 'Color', [0 0.4470 0.7410]); hold on;
    plot(R2Pm.Results{i}.Time, R2Pm.Results{i}.OutputIntensity, 'red*'); hold on;
    plot(R2Pm.Analysis{i}.Time,R2Pm.Analysis{i}.OutputIntensity,'k.','MarkerSize',20); hold on;
    title({'Resonant 2Photon Laser Power Measurement';  ['(', char(SystemProperties.FilePath.Data.Name(i)),'nm)']},'Interpreter','none')
    xlabel('Session Duration [min]'); ylabel('Laser Power Output Intensity [mW]')
    for j = 1:length(R2Pm.Analysis{i}.InterestValues)
         xline(R2Pm.Analysis{i}.Time(j),'-',{"Input Laser Power = " + num2str(R2Pm.Analysis{i}.InterestValues(j))+'%'})
    end
end

disp(R2Pm.Final.Final)
end