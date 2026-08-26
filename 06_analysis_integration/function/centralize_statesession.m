function [row, failure] = centralize_statesession(k, pax_fwhm, score, primary_analog, ...
    contents, stamp, old)
%CENTRALIZE_STATESESSION  One session's row of centralized_paxfwhm_state.
%   IN   k               1x1 double        which row of pax_table this is
%        pax_fwhm        1x1 struct        its centralized_paxfwhm data
%        score           1x1 struct | []   the sleep score, if it has one
%        primary_analog  1x1 struct | []   the whisker scoring instead
%        contents        Cx2 str           <property>, <field> the analysis runs over
%        stamp           1x1 double        what this row is rebuilt against
%        old             0|1 row table     what the last run wrote for this session
%   OUT  row             0x0 | 1x1 struct  the four fields are declared below
%        failure         1x1 str           "" unless the analysis threw
    row = struct('pax_i', {}, 'data', {}, 'input_stamp', {}, 'reused', {});
    failure = "";

    unmoved = ~isempty(old) && old.input_stamp == stamp;
    if ~unmoved
        try
            state_integrate = state_integration(score, primary_analog);
            state_integrate.trim_to_duration(pax_fwhm.t_axis(end));
            paxfwhm_state = state_linefwhm(state_integrate);
            paxfwhm_state.get_state_indices(pax_fwhm.t_axis, pax_fwhm.param.fs);
            for c = 1:size(contents, 1)
                property = contents(c, 1);
                field = contents(c, 2);
                if ~isfield(pax_fwhm, property) || ~isfield(pax_fwhm.(property), field)
                    continue
                end
                name = property + "_" + field;
                series = pax_fwhm.(property).(field);
                paxfwhm_state.get_summary(name, series);
                paxfwhm_state.get_powerdensity(name, series);
                paxfwhm_state.decompose_signal(name, series);
                paxfwhm_state.get_pppt_decomposition(name);
                paxfwhm_state.get_transitionsummary(name, series);
            end
            stripped = strip_state(paxfwhm_state, pax_fwhm.idx);
            row(1).pax_i = k;
            row(1).data = stripped;
            row(1).input_stamp = stamp;
            row(1).reused = false;
            return
        catch err
            % an analysis that throws is no reason to drop what the last run
            % computed, so the old row stands and the failure is named
            failure = string(err.message);
        end
    end
    if isempty(old)
        return
    end
    row(1).pax_i = k;
    row(1).data = old.data{1};
    row(1).input_stamp = stamp;
    row(1).reused = true;
end

function out = strip_state(paxfwhm_state, pax_idx)
    % The five result structs the downstream tables nest on, plus the two event
    % fields and the axes. sleep_obj is a handle to the state_integration that
    % built this, and carrying the object would carry the whole sleep score with
    % it; only the two pieces anything outside the class reads are kept --
    % rebuild_state_analysis takes param.transition_window off it, and the time
    % tables are what the windows mean.
    out = struct();
    for f = ["state_summary", "powerdensity", "band_decomposition", "transition", ...
            "peak_trough", "eventlist", "event", "t_axis", "state_idx", "param"]
        if isprop(paxfwhm_state, f)
            out.(f) = paxfwhm_state.(f);
        end
    end

    % The four boundary rows, NOT the nine series taken off them. bv, uppvs,
    % downpvs, totalpvs and eps are differences of these four and displacement is
    % each row against its own running median -- line_fwhm.m:133-137 and 168-186.
    % Keeping the source is half the size of keeping the derivatives, and it
    % cannot come to disagree with what line_fwhm computes from the same rows.
    out.idx = struct();
    for f = ["clean_upperBVboundary", "clean_lowerBVboundary", ...
            "clean_pvsupedge_idx", "clean_pvsdownedge_idx"]
        out.idx.(f) = pax_idx.(f);
    end

    out.sleep_obj = struct();
    if ~isempty(paxfwhm_state.sleep_obj)
        out.sleep_obj.param = paxfwhm_state.sleep_obj.param;
        out.sleep_obj.time_table = paxfwhm_state.sleep_obj.time_table;
        out.sleep_obj.info = paxfwhm_state.sleep_obj.info;
    end
end
