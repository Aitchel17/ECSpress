clc, clear
secondary_analysis_path = 'E:\OneDrive - The Pennsylvania State University\2023ecspress\02_secondary_analysis';
experiment_folder = 'G:\tmp\00_igkl';

% 1. Initialize TableManager
% Initialize Manager (Loads Excel automatically)
mtable_FWHM = tableManager(experiment_folder, secondary_analysis_path);
%% Table filtering
mtable_FWHM.filter_refTable("Primary_paxFWHM","paxfwhm");
%%
mtable_FWHM.filter_refTable("State_PaxFWHM","paxfwhm");
%%
mtable_FWHM.parseDepths();
%%
mtable_FWHM.filter_refTable("DepthLayer","L1");
%%
mtable_FWHM.filter_refTable("VesselID","PA");
%% 2. Aggregate Data to Vessels
% Loads FWHM data from .mat files and populates Vessel objects
mtable_FWHM.primaryTable = mtable_FWHM.aggregateData("Primary_paxFWHM");
mtable_FWHM.stateTable = mtable_FWHM.aggregateData("State_PaxFWHM");

%% Create Vessel Objects
Vessels = mtable_FWHM.aggregateVessels();
%%
v = Vessels(:,6)
%%
session_state = v.StateData(1);
transition_table = session_state.transition.thickness_bv;
%%
% Use the averaging function
averaged_transitions = average_transition_table(transition_table);
disp(averaged_transitions)


%% 3. Run Analysis
fprintf('Running average analysis on %d vessels...\n', numel(Vessels));
for v = 1:numel(Vessels)
    Vessels(v).computeAverage();
end
%%
% Select which transition state to plot (1 = first state, 2 = second state, etc.)
row_idx = 1;

figure()

% Collect mean_data(row_idx,:) from all vessels
all_vessel_traces_bv = [];
all_vessel_traces_pvs = [];
all_vessel_traces_eps = [];

for v_idx = 1:numel(Vessels)
    if ~isempty(Vessels(v_idx).AverageData) && isfield(Vessels(v_idx).AverageData, 'thickness_bv')
        avg_table_pvs = Vessels(v_idx).AverageData.thickness_totalpvs;
        avg_table_bv = Vessels(v_idx).AverageData.thickness_bv;
        avg_table_eps = Vessels(v_idx).AverageData.thickness_eps;
        % Check if row_idx exists in this vessel's table
        if height(avg_table_pvs) >= row_idx
            trace_bv = avg_table_bv.mean_data(row_idx,:);
            trace_pvs = avg_table_pvs.mean_data(row_idx,:);
            trace_eps = avg_table_eps.mean_data(row_idx,:);

            all_vessel_traces_bv = [all_vessel_traces_bv; trace_bv];
            all_vessel_traces_pvs = [all_vessel_traces_pvs; trace_pvs];
            all_vessel_traces_eps = [all_vessel_traces_eps; trace_eps];

        end
    end
end

% Get the state name from first vessel for title
state_name = Vessels(1).AverageData.thickness_totalpvs.state_name(row_idx);

% Compute grand average
grand_average_bv = mean(all_vessel_traces_bv, 1, 'omitnan');
grand_average_pvs = mean(all_vessel_traces_pvs, 1, 'omitnan');
grand_average_eps = mean(all_vessel_traces_eps, 1, 'omitnan');

% Plot
%plot(grand_average_bv, 'LineWidth', 2, Color = 'r')
hold on
%plot(grand_average_pvs, 'LineWidth', 2, Color = 'g')
 plot(grand_average_eps, 'LineWidth', 2)


% Optionally plot individual vessels in gray
for i = 1:size(all_vessel_traces_bv, 1)
    %plot(all_vessel_traces_bv(i,:), 'Color', [0.9 0.7 0.7])
    %plot(all_vessel_traces_pvs(i,:), 'Color', [0.7 0.9 0.7])
    plot(all_vessel_traces_eps(i,:), 'Color', [0.7 0.7 0.9])
end

xline(76, '--r', 'LineWidth', 1.5)
xlim([0 150])
xlabel('Data point')
ylabel('Thickness (pixel)')
title(sprintf('Average across %d Vessels - %s', size(all_vessel_traces_bv, 1), state_name))
% legend({'Grand Average', 'Individual Vessels'}, 'Location', 'best')
hold off
%%






%%
Vessels(1).AverageData.thickness_bv.mean_data(2,:)


%%
% Inspect results
if ~isempty(Vessels)
    disp('First Vessel Average Data:');
    disp(Vessels(1).AverageData);
end
