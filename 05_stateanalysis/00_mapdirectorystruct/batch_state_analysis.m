function batch_state_analysis(dirstruct_table, experiment_folder, primary_map, peripheral_map, stateanalysis_map, dirtable_dir, force_rerun)
%BATCH_STATE_ANALYSIS Automatically run state analysis for sessions that have
% paxfwhm.mat + sleep_score.mat but are missing paxfwhm_state.mat.
%
% Called from map_directory.m after write_dirtable.
%
% Inputs:
%   dirstruct_table   - Table from mapdirstruct (in-memory, already scanned)
%   experiment_folder - Root experiment folder (for re-scan at end)
%   primary_map       - primary_map cell array (for re-scan)
%   peripheral_map    - peripheral_map cell array (for re-scan)
%   stateanalysis_map - stateanalysis_map cell array (for re-scan)
%   dirtable_dir      - Path to dirtable .xlsx (for re-scan update)
%   force_rerun       - (optional) true = 기존 paxfwhm_state.mat도 강제 재처리 (default: false)

if nargin < 7 || isempty(force_rerun)
    force_rerun = false;
end

% --- Filter: paxfwhm + sleep_score 있고, state analysis 없는 행 ---
has_paxfwhm    = string(dirstruct_table.Primary_paxFWHM)        == "paxfwhm.mat";
has_sleepscore = string(dirstruct_table.Peripheral_SleepScore)  == "sleep_score.mat";
no_state       = string(dirstruct_table.State_PaxFWHM)          == "NA";

if force_rerun
    target_rows = find(has_paxfwhm & has_sleepscore);  % 기존 것도 포함
    fprintf('[batch_state_analysis] force_rerun=true: %d session(s) queued (including existing).\n', numel(target_rows));
else
    target_rows = find(has_paxfwhm & has_sleepscore & no_state);  % 신규만
    fprintf('[batch_state_analysis] %d session(s) queued (new only).\n', numel(target_rows));
end

if isempty(target_rows)
    fprintf('[batch_state_analysis] Nothing to process. All sessions are up to date.\n');
    return;
end

% Contents to process (matches stateintegration_main.m)
contents_types = {["thickness","eps","bv","totalpvs","dynamic_pvs","static_pvs"], ...
    ["displacement","dynamicpvs","staticpvs","dynamicbv","staticbv"]};

for row_i = 1:numel(target_rows)
    sessiondir = char(dirstruct_table.Directory{target_rows(row_i)});
    fprintf('[%d/%d] %s\n', row_i, numel(target_rows), sessiondir);

    try
        %% 1. Load paxfwhm via ECSSession (matches stateintegration_main.m)
        % ECSSession으로 로드해야 pax_fwhm 객체 구조가 올바르게 복원됨
        session = ECSSession(sessiondir);
        session = session.load_primary_results();
        pax_fwhm = session.pax_fwhm;

        %% 2. Build state_integration (auto-finds sleep_score.mat in peripheral/)
        state_integrate = state_integration(sessiondir);

        % 이미징 중단으로 pax_fwhm이 sleep_score보다 짧을 수 있음
        % -> 짧은 쪽(이미징 끝시간)에 맞춰 time_table 클립
        state_integrate.trim_to_duration(pax_fwhm.t_axis(end));

        %% 3. Build state_linefwhm and attach time axis
        paxfwhm_state = state_linefwhm(state_integrate);
        paxfwhm_state.get_state_indices(pax_fwhm.t_axis, pax_fwhm.param.fs);

        %% 4. Run analysis for each content type
        % isfield()은 struct에만 동작 -> class 객체에는 try/catch로 안전하게 접근
        for typeidx = 1:numel(contents_types)
            contents_type = contents_types{typeidx};
            d_type = contents_type(1);
            for cidx = 2:numel(contents_type)
                content_name = contents_type(cidx);
                name = strcat(d_type, '_', content_name);

                % try/catch: struct이든 class 객체든 필드 없으면 스킵
                try
                    data = pax_fwhm.(char(d_type)).(char(content_name));
                catch
                    fprintf('    [skip] %s (field not found)\n', name);
                    continue;
                end

                fprintf('    %s\n', name);
                paxfwhm_state.get_summary(name, data);
                paxfwhm_state.get_powerdensity(name, data);
                paxfwhm_state.decompose_signal(name, data);
                paxfwhm_state.get_pppt_decomposition(name);
                paxfwhm_state.get_transitionsummary(name, data);
            end
        end

        %% 5. Save paxfwhm_state.mat -> session/state_analysis/
        paxfwhm_state.save2disk('paxfwhm_state', state_integrate.dir_struct.stateanalysis);
        fprintf('  -> Saved: %s\n', fullfile(state_integrate.dir_struct.stateanalysis, 'paxfwhm_state.mat'));

    catch ME
        warning('  [FAILED] %s\n  Reason: %s', sessiondir, ME.message);
    end
end

%% 6. Re-scan dirtable to reflect newly created paxfwhm_state.mat files
if nargin >= 6 && ~isempty(dirtable_dir)
    fprintf('[batch_state_analysis] Re-scanning dirtable...\n');
    updated_table = mapdirstruct(experiment_folder, primary_map, peripheral_map, stateanalysis_map);
    write_dirtable(updated_table, dirtable_dir);
    fprintf('[batch_state_analysis] dirtable updated: %s\n', dirtable_dir);
end

end
