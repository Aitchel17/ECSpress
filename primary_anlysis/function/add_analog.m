function primary_analog = add_analog(loaded_data,mdfExtract)
%ADD_ANALOG Summary of this function goes here
%   Detailed explanation goes here
% load and process analog channel if analog data not exist in loaded data
if ~isfield(loaded_data, 'analog')
    analog = mdfExtract.loadanalog;
    primary_analog = analysis_analog(analog.info,analog.data); % analysis analog class construct
    if ~isfield(primary_analog, 'raw_Air_puff1') && isfield(analog.data, 'raw_Air_puff1')
        primary_analog.airtable = primary_analog.get_airtable('raw_Air_puff1');
    end
    if ~isfield(primary_analog, 'ecogspectrum') && isfield(analog.data, 'raw_ECoG')
        primary_analog.ecogspectrum = primary_analog.get_ecogspectrum('raw_ECoG');
    end
else
    disp('loading previously saved analog data')
    primary_analog = loaded_data.analog;
    % ensure
    if ~isfield(primary_analog, 'airtable') && isfield(primary_analog.data, 'raw_Air_puff1')
        disp('Get airtable')
        primary_analog.airtable = primary_analog.get_airtable('raw_Air_puff1');
    end
    if ~isfield(primary_analog, 'ecogspectrum') && isfield(primary_analog.data, 'raw_ECoG')
        disp('Calculate ecogspectrum')
        primary_analog.ecogspectrum = primary_analog.get_ecogspectrum('raw_ECoG');
    end
end
end

