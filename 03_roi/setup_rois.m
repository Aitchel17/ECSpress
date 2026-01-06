function setup_rois(roilist, twophoton_processed)
% SETUP_ROIS Defines the standard set of manual ROIs for the experiment.
%   Establishes a hierarchy of manual ROIs with copy/modify logic.
%   Prompts user for channel selection when modifying.

%% 1. Manual Extraparenchyma (Base) -> Ch2
process_roi_step(roilist, twophoton_processed, 'manual_extraparenchyma', '', 'polygon', 2);

%% 2. Dilated PVS (Copy from Extraparenchyma) -> Ch2
process_roi_step(roilist, twophoton_processed, 'manual_dilated_pvs', 'manual_extraparenchyma', 'polygon', 2);

%% 3. Constricted PVS (Copy from Dilated PVS) -> Ch2
process_roi_step(roilist, twophoton_processed, 'manual_constricted_pvs', 'manual_dilated_pvs', 'polygon', 2);

%% 4. Dilated BV (Copy from Dilated PVS) -> Ch1
process_roi_step(roilist, twophoton_processed, 'manual_dilated_bv', 'manual_dilated_pvs', 'polygon', 1);

%% 5. Constricted BV (Copy from Constricted PVS) -> Ch1
process_roi_step(roilist, twophoton_processed, 'manual_constricted_bv', 'manual_constricted_pvs', 'polygon', 1);

%% 6. Lines
process_roi_step(roilist, twophoton_processed, 'manual_dipin', '', 'line', 1);
process_roi_step(roilist, twophoton_processed, 'manual_dipout', '', 'line', 2);

% Save
roilist.save2disk();

end

function process_roi_step(roilist, data, target_name, source_name, mode, default_ch_idx)
% Helper to handle check -> copy -> modify logic

current_labels = roilist.list();
exists = any(strcmp(current_labels, target_name));

if exists
    % ROI Exists: Ask to modify
    fprintf('\nROI "%s" already exists.\n', target_name);
    prompt = sprintf('Do you want to MODIFY "%s"? [y/N]: ', target_name);
    user_choice = input(prompt, 's');

    if strcmpi(user_choice, 'y')
        % Ask for channel
        fprintf('Select channel to display for modification:\n');
        fprintf('1: Channel 1 (BV)\n');
        fprintf('2: Channel 2 (PVS)\n');
        ch_choice = input(sprintf('Choice [Default %d]: ', default_ch_idx));

        if isempty(ch_choice)
            ch_choice = default_ch_idx;
        end

        if ch_choice == 1
            img = data.ch1;
        else
            img = data.ch2; % Default to Ch2 if invalid or 2
        end

        roilist.modifyroi(img, target_name);
    end

else
    % ROI Missing: Create (Copy or New)
    fprintf('\nCreating NEW ROI: "%s"...\n', target_name);

    copied = false;

    % Try to copy from source if provided
    if ~isempty(source_name)
        if any(strcmp(current_labels, source_name))
            try
                roilist.copyroi(source_name, target_name);
                fprintf('  -> Copied successfully from "%s".\n', source_name);
                copied = true;
            catch ME
                fprintf('  -> Failed to copy from "%s": %s\n', source_name, ME.message);
            end
        else
            fprintf('  -> Source "%s" not found. Cannot copy.\n', source_name);
        end
    end

    % Select Image for creation/modification defaults
    if default_ch_idx == 1
        img = data.ch1;
    else
        img = data.ch2;
    end

    if copied
        % If copied, enter modify mode immediately so user can adjust it
        fprintf('  -> Opening for adjustment...\n');
        roilist.modifyroi(img, target_name);
    else
        % If not copied, create from scratch
        fprintf('  -> Drawing from scratch...\n');
        roilist.addormodifyroi(img, target_name, mode);
    end

end
end
