function results = get_all_spectra(bouts, analysis_obj, fs)
%GET_ALL_SPECTRA Calculate spectrograms for all 1D data in analysis object
%   bouts: cell array of indices
%   analysis_obj: e.g. session.pax_fwhm (must have t_axis)
%   fs: sampling frequency

results = struct();

if ~isprop(analysis_obj, 't_axis') && ~isfield(analysis_obj, 't_axis')
    warning('Analysis object missing t_axis. Skipping spectral analysis.');
    return;
end

t_axis_len = numel(analysis_obj.t_axis);
props = fieldnames(analysis_obj);
if isobject(analysis_obj)
    props = properties(analysis_obj);
end

for i = 1:numel(props)
    propName = props{i};
    if strcmp(propName, 't_axis') || strcmp(propName, 'param'), continue; end

    targetStruct = analysis_obj.(propName);

    if isstruct(targetStruct)
        subFields = fieldnames(targetStruct);

        for k = 1:numel(subFields)
            fName = subFields{k};
            data = targetStruct.(fName);

            % Check if data is 1D and matches time length
            if isnumeric(data) && isvector(data) && (numel(data) == t_axis_len)
                % 1. Compute Spectrogram List
                spect_list = get_spectboutarray(bouts, data, fs, 0.1);

                % 2. Compute Summary
                summary = get_summaryspectrogram(spect_list, fs);

                % 3. Decompose Signal (Continuous, VLF, LF)
                try
                    decomposed_full = decompose_signal(data, fs);
                    band_names = fieldnames(decomposed_full);

                    decomposed_bouts = struct();
                    for b_idx = 1:numel(band_names)
                        bn = band_names{b_idx};
                        full_band_signal = decomposed_full.(bn);

                        % Slice by bouts
                        cell_bouts = cell(size(bouts));
                        for c = 1:numel(bouts)
                            cell_bouts{c} = full_band_signal(bouts{c});
                        end
                        decomposed_bouts.(bn) = cell_bouts;
                    end
                catch ME
                    warning('Decomposition failed for %s.%s: %s', propName, fName, ME.message);
                    decomposed_bouts = [];
                end

                % Store
                results.(propName).(fName).list = spect_list;
                results.(propName).(fName).summary = summary;
                results.(propName).(fName).decomposed = decomposed_bouts;
            end
        end
    end
end
end
