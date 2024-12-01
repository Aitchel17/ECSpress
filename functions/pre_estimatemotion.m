function [pixelShift_table] = pre_estimatemotion(dft_stack,ref_frame,Vertices)
    disp(['Estimate motion to get drift table using reference frame ' num2str(ref_frame)]);
    xy = round(Vertices);

    % Extract the selected region from the stack
    hold_stack = dft_stack(xy(1,2):xy(3,2), xy(1,1):xy(3,1), :);
    
    % Display the selected region in a slice viewer
    figure(5);
    sliceViewer(hold_stack);
    
    % Perform Fourier Transform on the first slice (used as reference)
    first_fft = fft2(hold_stack(:,:,ref_frame));
    
    % Initialize the table to store pixel shifts
    pixelShift_table = zeros(4, size(hold_stack, 3));  % 4 rows for shift values (x, y, and shifts)
    
    % Loop over all slices in the stack
    for sli = 1:size(hold_stack, 3)
        
        
        % Fourier transform of the current slice
        regframe = fft2(hold_stack(:,:,sli));
        
        % Estimate the pixel shift using DFT registration
        [pixelShift_table(:, sli), ~] = dft_registration(first_fft, regframe);
    end
    
    % Display the X and Y shifts
    figure('name', 'Pixel Shift', 'NumberTitle', 'off');
    subplot(2, 1, 1);
    xshift = medfilt1(pixelShift_table(4,:), 100);  % Apply median filter to smooth x-shifts
    plot(xshift);
    title('X shift');
    
    subplot(2, 1, 2);
    yshift = medfilt1(pixelShift_table(3,:), 100);  % Apply median filter to smooth y-shifts
    plot(yshift);
    title('Y shift');
end

