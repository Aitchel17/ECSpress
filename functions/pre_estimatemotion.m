function [pixelShift_table] = pre_estimatemotion(dst_stack,ref_frame,Vertices)
    raw_hold_image=double(dst_stack);
    hold_stack=raw_hold_image(Vertices(1,2):Vertices(3,2), Vertices(1,1):Vertices(3,1),:);
    first_fft = fft2(ref_frame);

    pixelShift_table = [];
    for sli = drange(1:size(hold_stack,3))
        disp(sli)
        regframe = fft2(hold_stack(:,:,sli));
        [pixelShift_table(:,sli),~] = dft_registration(first_fft,regframe);
% Inputs
% buf1ft    Fourier transform of reference image, 
%           DC in (1,1)   [DO NOT FFTSHIFT]
% buf2ft    Fourier transform of image to register, 
%           DC in (1,1) [DO NOT FFTSHIFT]
% usfac     Upsampling factor (integer). Images will be registered to 
%           within 1/usfac of a pixel. For example usfac = 20 means the
%           images will be registered within 1/20 of a pixel. (default = 1)
%
% Outputs
% output =  [error,diffphase,net_y_shift,net_x_shift]
% error     Translation invariant normalized RMS error between f and g
% diffphase     Global phase difference between the two images (should be
%               zero if images are non-negative).
% net_row_shift net_col_shift   Pixel shifts between images
% Greg      (Optional) Fourier transform of registered version of buf2ft,
%           the global phase difference is compensated for.
    end

    figure('name', 'pixel shift','NumberTitle','off')
    subplot(2,1,1)
    xshift = pixelShift_table(4,:);
    plot(xshift)
    title('X shift')
    subplot(2,1,2)
    yshift = pixelShift_table(3,:);
    plot(yshift)
    title('Y-shift')
end

