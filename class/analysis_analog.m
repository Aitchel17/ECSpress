classdef analysis_analog
    %ANALYSIS_ANALOG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        info = struct()
        data = struct()
        airtable = table()
        ecogspectrum = struct()
    end
    
    methods
        function obj = analysis_analog(info_struct,data_struct)
            %ANALYSIS_ANALOG Construct an instance of this class
            %   Detailed explanation goes here
            obj.info = info_struct;
            obj.data = data_struct;
            obj.data.taxis = linspace(0,str2double(obj.info.analogcount)/str2double(obj.info.analogfreq(1:end-3)),str2double(obj.info.analogcount));
        end
        
        function air_puff_table = get_airtable(obj,airpuff_fieldname)
            % find single air puff initiation
            binary_airpuff = obj.data.(airpuff_fieldname)>0; % binarize data on or off
            diff_airpuff = diff(binary_airpuff); % differentiate rising and faling
            stim_on_idx = find(diff_airpuff ==1); % find rising edge, each data point is cumulative point when stim on
            stim_off_idx = find(diff_airpuff ==-1); % find falling edge
            stim_on_time = obj.data.taxis(stim_on_idx); % match time scale, each data point is time when stim on
            stim_off_time = obj.data.taxis(stim_off_idx+1);
            stim_on_int = diff(stim_on_time); % slope of stim on time point = frequency of stim
            % find session boundary, inter session interval should be
            % longer then 10 sec
            session_boundary_idx = find(stim_on_int>10); % the interval between session should longer then 10 sec
            session_end = [stim_off_time(session_boundary_idx),stim_off_time(end)];
            session_start = [stim_on_time(1),stim_on_time(session_boundary_idx+1)];
            session_duration = session_end-session_start;
            stim_on_off_idx = stim_off_idx-stim_on_idx;
            stim_on_on_idx = stim_on_idx(2:end)-stim_on_idx(1:end-1);
            stim_duty = stim_on_off_idx(1:end-1)./stim_on_on_idx;
            session_frequency = [];
            session_duty = [];
            for session_id = 0:length(session_boundary_idx)
            % initial
                if session_id == 0
                    session_start_idx = 1;
                    session_end_idx = session_boundary_idx(1)-1;
                elseif session_id == length(session_boundary_idx)
                    session_start_idx = session_boundary_idx(session_id)+1;
                    session_end_idx = length(stim_on_int);
                else
                    session_start_idx = session_boundary_idx(session_id)+1;
                    session_end_idx = session_boundary_idx(session_id+1)-1;
                end
                session_frequency = [session_frequency,1/mean(stim_on_int(session_start_idx:session_end_idx))];
                session_duty = [session_duty,mean(stim_duty(session_start_idx:session_end_idx))];
            end
            air_puff_data = [session_start;session_end;session_duration;session_frequency;session_duty]';
            air_puff_table = array2table(air_puff_data,'VariableNames', {'StartTime', 'EndTime', 'Duration', 'FrequencyHz', 'DutyCycle'});
        end

        function ecog_spectrum = get_ecogspectrum(obj,ecog_fieldname, samplingfrequency)
            arguments
                obj
                ecog_fieldname
                samplingfrequency = str2double(obj.info.analogfreq(1:end-2));
            end
            % ECoG processing
            ECoG = obj.data.(ecog_fieldname);
            [ecog_spectrum] = analog_ecogspectrum(samplingfrequency,ECoG);
        
        end
    end

end

