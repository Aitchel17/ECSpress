function bout_struct = add_state_analysis(analysis, varargin)
%ADD_STATE_ANALYSIS Slice analysis data into sleep state bouts
%   Input: analysis object (must have t_axis), timebins (e.g., NREMTimes)

bout_struct = struct();

if strcmp(class(analysis), 'line_fwhm')
    if ~isprop(analysis, 't_axis') && ~isfield(analysis, 't_axis')
        error('Analysis object must have t_axis property or field');
    end
    t_axis = analysis.t_axis;
    data_len = numel(t_axis);
    property_fieldnames = fieldnames(analysis);

    % 1. Calculate Indices and Bouts
    stateidx = statebin2frame(varargin{:}, t_axis);
    bouts = stateidx2bouts(stateidx);

    % 2. Initialize Output
    bout_struct.t_axis = t_axis;
    bout_struct.stateidx = stateidx;
    bout_struct.bouts = bouts;
    if isprop(analysis, 'param') || isfield(analysis, 'param')
        bout_struct.param = analysis.param;
        for property_idx = 1:numel(property_fieldnames)
            for cnt  = 1:numel(bouts)
                data_len = numel(analysis.t_axis);

                target_property = analysis.(property_fieldnames{property_idx});
                if isstruct(target_property)
                    entry = struct();
                    entry.bout_idx = cnt;
                    contents_names = fieldnames(target_property);
                    for idx = 1:numel(contents_names)
                        data = analysis.(property_fieldnames{property_idx}).(contents_names{idx});

                        if isnumeric(data) || islogical(data)
                            sz = size(data);
                            t_dim = find(sz == data_len, 1, 'last'); % Assume time is last matching dim if ambiguous

                            if ~isempty(t_dim)
                                % Prepare dynamic slicing indices
                                idx_subs = repmat({':'}, 1, ndims(data));
                                idx_subs{t_dim} = bouts{cnt};
                                entry.(contents_names{idx}) = data(idx_subs{:});
                            end
                        end
                    end
                    if ~strcmp(property_fieldnames{property_idx},'param') && (numel(fieldnames(entry)) > 1)
                        bout_struct.(property_fieldnames{property_idx})(cnt) = entry;
                    end
                end
            end
        end
    end
end
end