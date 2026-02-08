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
        subFields_names = fieldnames(targetStruct);

        for k = 1:numel(subFields_names)
            fName = subFields_names{k};
            data = targetStruct.(fName);

            % Check if data is 1D and matches time length
            if isnumeric(data) && isvector(data) && (numel(data) == t_axis_len)
                % 1. Compute Spectrogram List
                spect_list = get_spectboutarray(bouts, data, fs, 0.1);

                % 2. Compute Summary
                spec_summary = get_summaryspectrogram(spect_list, fs);

                % 3. Decompose Signal (Continuous, VLF, LF)
                try
                    decomposed_full = decompose_signal(data, fs);
                    band_names = fieldnames(decomposed_full);

                    decomposed_bouts = struct();
                    peak_analysis = struct();
                    for b_idx = 1:numel(band_names)
                        bandname = band_names{b_idx};
                        full_band_signal = decomposed_full.(bandname);

                        % Slice by bouts
                        % Slice by bouts
                        n_bouts = numel(bouts);
                        cell_bouts = cell(n_bouts, 1);

                        % Vectors -> Cell arrays
                        p2p_time_cell = cell(n_bouts, 1);
                        p2t_amp_cell = cell(n_bouts, 1);
                        peak_idx_cell = cell(n_bouts, 1);
                        peak_val_cell = cell(n_bouts, 1);
                        trough_idx_cell = cell(n_bouts, 1);
                        trough_val_cell = cell(n_bouts, 1);

                        % Scalars -> Arrays
                        p2p_avg_arr = nan(n_bouts, 1);
                        p2p_std_arr = nan(n_bouts, 1);
                        p2t_avg_arr = nan(n_bouts, 1);
                        p2t_std_arr = nan(n_bouts, 1);
                        median_arr =  nan(n_bouts, 1);
                        median_peak_arr = nan(n_bouts, 1);
                        median_trough_arr = nan(n_bouts, 1);


                        for bout_idx = 1:numel(bouts)
                            cell_bouts{bout_idx} = full_band_signal(bouts{bout_idx});
                            % Peak Analysis
                            p_stats = get_peak_analysis(cell_bouts{bout_idx}, fs);

                            % Store Vectors/Raw Data
                            p2p_time_cell{bout_idx} = p_stats.peak2peak_time;
                            p2t_amp_cell{bout_idx} = p_stats.peak2trough_amp;
                            peak_idx_cell{bout_idx} = p_stats.peak_idx;
                            peak_val_cell{bout_idx} = p_stats.peak_val;
                            trough_idx_cell{bout_idx} = p_stats.trough_idx;
                            trough_val_cell{bout_idx} = p_stats.trough_val;
                            median_arr(bout_idx) = median(cell_bouts{bout_idx});
                            median_peak_arr = median(p_stats.peak_val);
                            median_trough_arr = median(p_stats.trough_val);
                            
                            % Store Scalars
                            p2p_avg_arr(bout_idx) = p_stats.p2p_avg;
                            p2p_std_arr(bout_idx) = p_stats.p2p_std;
                            p2t_avg_arr(bout_idx) = p_stats.p2t_avg;
                            p2t_std_arr(bout_idx) = p_stats.p2t_std;
                        end
                        decomposed_bouts.(bandname).data = cell_bouts;
                        decomposed_bouts.(bandname).median = median_arr(bout_idx);

                        % Assign to structure (Vectors)
                        peak_analysis.(bandname).peak2peak_time = p2p_time_cell;
                        peak_analysis.(bandname).peak2trough_amp = p2t_amp_cell;
                        peak_analysis.(bandname).peak_idx = peak_idx_cell;
                        peak_analysis.(bandname).peak_val = peak_val_cell;
                        peak_analysis.(bandname).trough_idx = trough_idx_cell;
                        peak_analysis.(bandname).trough_val = trough_val_cell;

                        % Assign to structure (Scalars)
                        peak_analysis.(bandname).p2p_avg = p2p_avg_arr;
                        peak_analysis.(bandname).p2p_std = p2p_std_arr;
                        peak_analysis.(bandname).p2t_avg = p2t_avg_arr;
                        peak_analysis.(bandname).p2t_std = p2t_std_arr;
                        peak_analysis.(bandname).peak_median = median_peak_arr;
                        peak_analysis.(bandname).trogh_median = median_trough_arr;
                    end
                catch ME
                    warning('Decomposition failed for %s.%s: %s', propName, fName, ME.message);
                    decomposed_bouts = [];
                    peak_analysis = [];
                end

                % Store
                results.(propName).(fName).list = spect_list;
                results.(propName).(fName).spec_summary = spec_summary;
                results.(propName).(fName).decomposed = decomposed_bouts;
                results.(propName).(fName).peak_analysis = peak_analysis;
            end
        end
    end
end
end
