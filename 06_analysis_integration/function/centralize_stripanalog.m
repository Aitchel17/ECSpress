function out = centralize_stripanalog(primary_analog)
%CENTRALIZE_STRIPANALOG  What centralized_analysis_analog keeps out of an analysis_analog.
%   IN   primary_analog  1x1 analysis_analog
%   OUT  out             1x1 struct

    % The resampled traces and the normalised spectrogram; nothing at the original
    % sample rate. rawdata, emg.power/signal and ball's pre-resample traces are all
    % 1.8 M samples per field and nothing downstream reads them. Of the four arrays
    % in ecogspectrum only baseline_S is taken -- it is the log ratio that is
    % actually looked at, and S, log_norm_spectrum and awake_S are the same
    % measurement in three other framings.
    %
    % Coming out as a plain struct also drops the dependency on the analysis_analog
    % classdef, which lives in 02_othersignal rather than either repo root here.
    % Fields are guarded because a sleep recording has no airtable or ball and a
    % whisker one has no sleep scoring behind its spectrogram.
    out = struct();
    if isprop(primary_analog, 'airtable')
        out.airtable = primary_analog.airtable;
    end
    out.ball = pick_fields(primary_analog, 'ball', ...
        ["ds_fps", "rs_taxis", "resampled_ball"]);
    out.ecog = pick_fields(primary_analog, 'ecog', ...
        ["ds_fps", "resampled_ecog", "resampled_taxis"]);
    out.emg = pick_fields(primary_analog, 'emg', ...
        ["ds_fps", "rs_taxis", "resampled_power", "resampled_signal"]);

    out.ecogspectrum = struct();
    if isprop(primary_analog, 'ecog') && isfield(primary_analog.ecog, 'ecogspectrum')
        spectrum = primary_analog.ecog.ecogspectrum;
        for f = ["T", "F", "params", "movingwin"]
            if isfield(spectrum, f)
                out.ecogspectrum.(f) = spectrum.(f);
            end
        end
        % baseline_S is the log ratio against an awake baseline, so it only exists
        % where there was sleep scoring to define that baseline. A whisker
        % recording has none and keeps the raw power S instead. The two are NOT
        % the same quantity and are not even the same way round -- baseline_S is
        % frequency by time, S is time by frequency -- so they keep their own
        % names and a reader has to look at which one is there.
        if isfield(spectrum, 'baseline_S')
            out.ecogspectrum.baseline_S = single(spectrum.baseline_S);
        elseif isfield(spectrum, 'S')
            out.ecogspectrum.S = single(spectrum.S);
        end
    end
end

function out = pick_fields(container, group, keep)
    out = struct();
    if ~isprop(container, group) || isempty(container.(group))
        return
    end
    for f = keep
        if isfield(container.(group), f)
            out.(f) = container.(group).(f);
        end
    end
end
