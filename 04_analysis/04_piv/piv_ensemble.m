function [xtable, ytable, utable, vtable, typevector, correlation_map, result_conv_ensemble] = piv_ensemble(imgstack, window_sizes, subpixfinder, roi_inpt, mask_inpt, imdeform, repeat, mask_auto, do_pad, use_gpu)
%PIV_ENSEMBLE  Ensemble multipass FFT PIV on an in-memory image stack.
%   Logic is equivalent to piv_FFTensemble. Deviations are marked [DEVIATION].
%
%   INPUTS
%     imgstack     : double H x W x N array (pre-processed, values in [0,1])
%     window_sizes : P x 2 matrix [window_size, step] per pass, or P x 1.
%                    e.g. [64 32; 32 16; 16 8] for 3 passes. Max 4 rows.
%     subpixfinder : 1 = 3-point Gaussian, 2 = 2D Gaussian
%     roi_inpt     : [x y w h] ROI, or [] for full frame
%     mask_inpt    : H x W mask (1=masked), or []
%     imdeform     : image deformation method, e.g. '*spline'
%     repeat       : 1 = repeated/multishift correlation on last pass
%     mask_auto    : 1 = suppress auto-correlation peak
%     do_pad       : 1 = linear zero-padded correlation on last pass
%
%   OUTPUTS
%     xtable, ytable       : vector grid coordinates [pixels, global frame]
%     utable, vtable       : displacement [pixels/frame]
%     typevector           : 1=valid, 0=masked
%     correlation_map      : normalised correlation strength
%     result_conv_ensemble : accumulated correlation volume [ia x ia x Nvec]

warning off %#ok<WNOFF>

%% [DEVIATION] Parse consolidated window_sizes argument.
% piv_FFTensemble receives: interrogationarea, step, passes, int2, int3, int4 separately.
if size(window_sizes,1) > 4
    warning('piv_ensemble: max 4 passes. Truncating.');
    window_sizes = window_sizes(1:4,:);
end
passes = size(window_sizes,1);
ws = window_sizes(:,1);
if size(window_sizes,2) >= 2; ss = window_sizes(:,2); else; ss = floor(ws/2); end
interrogationarea = ws(1);
step              = ss(1);
int2 = ws(min(2,end)); int3 = ws(min(3,end)); int4 = ws(min(4,end));

%% [DEVIATION] Input is an in-memory H x W x N matrix, not file paths / VideoReader.
% piv_FFTensemble: amount_input_imgs = size(filepath,1) or numel(video_frame_selection)
% piv_FFTensemble processes independent pairs (1-2, 3-4, ...).
% Here we process all consecutive pairs (1-2, 2-3, ...) for richer ensemble averaging.
[H, W, N] = size(imgstack);
n_pairs = N - 1;
if n_pairs < 1; error('piv_ensemble: need at least 2 frames.'); end

if use_gpu
    imgstack = gpuArray(single(imgstack));
    mask_inpt = gpuArray(mask_inpt);
    predictor_interp_method = 'cubic';
    linear_interp_method = 'linear';
    if contains(imdeform, 'spline', 'IgnoreCase', true)
        imdeform = 'cubic';
    else
        imdeform = strrep(imdeform, '*', '');
    end
else
    predictor_interp_method = '*spline';
    linear_interp_method = '*linear';
end

%% [DEVIATION] No GUI, no cancel, no autolimit, no bg subtraction, no PIVlab_preproc.
% Preprocessing is done externally via opticalflow_preprocess.


%% [DEVIATION] Mask setup — piv_FFTensemble builds average_mask from converted_mask
% cell array (one entry per pair). Here we use a single static mask_inpt and
% accumulate average_mask_roi incrementally across the frame loop.
if isempty(mask_inpt)
    mask_inpt = zeros(H, W); 
end


%% Pass 1: ensemble accumulation — piv_FFTensemble lines 7-357
result_conv_ensemble = zeros(interrogationarea, interrogationarea); % piv_FFTensemble line 7
corr_map_cnt = 0;
%% ROI setup — piv_FFTensemble lines 82-98
if numel(roi_inpt) > 0
    xroi=roi_inpt(1);
    yroi=roi_inpt(2);
    widthroi=roi_inpt(3);
    heightroi=roi_inpt(4);
else
    xroi=0; 
    yroi=0; 
    widthroi=W-1; 
    heightroi=H-1;
end
average_mask_roi = zeros(heightroi+1, widthroi+1);


for k = 1:n_pairs % [DEVIATION] piv_FFTensemble: for ensemble_i1=1:2:amount_input_imgs
    if mod(k, 10) == 0; fprintf('.'); if use_gpu; wait(gpuDevice); end; end

    %% [DEVIATION] Load images from imgstack instead of file / VideoReader.
    % piv_FFTensemble also applies: RGB→grey, bg subtraction, stretchlim, PIVlab_preproc.
    image1_roi = double(imgstack(yroi:yroi+heightroi, xroi:xroi+widthroi, k));
    image2_roi = double(imgstack(yroi:yroi+heightroi, xroi:xroi+widthroi, k+1));
    mask_inpt_roi    = mask_inpt(yroi:yroi+heightroi, xroi:xroi+widthroi);
    average_mask_roi = average_mask_roi + mask_inpt_roi; % [DEVIATION] incremental
    
    % piv_FFTensemble lines 99-142

    mask = mask_inpt_roi; 
    gen_image1_roi = image1_roi; 
    gen_image2_roi = image2_roi;

    %% Grid on UNPADDED image
    miniy=1+(ceil(interrogationarea/2));
    minix=1+(ceil(interrogationarea/2));
    maxiy=step*(floor(size(image1_roi,1)/step))-(interrogationarea-1)+(ceil(interrogationarea/2));
    maxix=step*(floor(size(image1_roi,2)/step))-(interrogationarea-1)+(ceil(interrogationarea/2));

    numelementsy=floor((maxiy-miniy)/step+1);
    numelementsx=floor((maxix-minix)/step+1);

    LAy=miniy; LAx=minix;
    LUy=size(image1_roi,1)-maxiy; LUx=size(image1_roi,2)-maxix;
    shift4centery=round((LUy-LAy)/2); shift4centerx=round((LUx-LAx)/2);
    if shift4centery<0; shift4centery=0; end
    if shift4centerx<0; shift4centerx=0; end
    miniy=miniy+shift4centery; minix=minix+shift4centerx;
    maxix=maxix+shift4centerx; maxiy=maxiy+shift4centery;
    %% Pad AFTER grid computation
    image1_roi=padarray(image1_roi,[ceil(interrogationarea/2) ceil(interrogationarea/2)],min(min(image1_roi)));
    image2_roi=padarray(image2_roi,[ceil(interrogationarea/2) ceil(interrogationarea/2)],min(min(image1_roi)));
    mask      =padarray(mask,      [ceil(interrogationarea/2) ceil(interrogationarea/2)],0);
    if (rem(interrogationarea,2)==0);  % SubPixOffset
        SubPixOffset=1; 
    else
        SubPixOffset=0.5
    end

    if k==1     % Initialise tables on first pair
        xtable=zeros(numelementsy,numelementsx); ytable=xtable;
        utable=xtable; vtable=xtable;
        typevector=ones(numelementsy,numelementsx);
        correlation_map=zeros(size(typevector));
    end

    %% Sub-image index arrays — piv_FFTensemble lines 169-176
    s0=(repmat((miniy:step:maxiy)'-1,1,numelementsx)+repmat(((minix:step:maxix)-1)*size(image1_roi,1),numelementsy,1))';
    s0=permute(s0(:),[2 3 1]);
    s1=repmat((1:interrogationarea)',1,interrogationarea)+repmat(((1:interrogationarea)-1)*size(image1_roi,1),interrogationarea,1);
    ss1=repmat(s1,[1,1,size(s0,3)])+repmat(s0,[interrogationarea,interrogationarea,1]);

    image1_cut=image1_roi(ss1);
    image2_cut=image2_roi(ss1);

    %% do_pad (pass 1 only) — piv_FFTensemble lines 
    if do_pad==1 && passes==1
        image1_cut=image1_cut-mean(image1_cut,[1 2]);
        image2_cut=image2_cut-mean(image2_cut,[1 2]);
        image1_cut=[image1_cut zeros(interrogationarea,interrogationarea-1,size(image1_cut,3));zeros(interrogationarea-1,2*interrogationarea-1,size(image1_cut,3))];
        image2_cut=[image2_cut zeros(interrogationarea,interrogationarea-1,size(image2_cut,3));zeros(interrogationarea-1,2*interrogationarea-1,size(image2_cut,3))];
    end

    %% FFT correlation — piv_FFTensemble line 179
    result_conv=fftshift(fftshift(real(ifft2(conj(fft2(image1_cut)).*fft2(image2_cut))),1),2);
    if do_pad==1 && passes==1
        result_conv=result_conv((interrogationarea/2):(3*interrogationarea/2)-1,(interrogationarea/2):(3*interrogationarea/2)-1,:);
    end

	%% repeated  Correlation in the first pass (might make sense to repeat more often to make it even more robust...)
	if repeat == 1 && passes == 1
		ms=round(step/4); %multishift parameter so groß wie viertel int window
		%Shift left bot
		s0B = (repmat((miniy+ms:step:maxiy+ms)'-1, 1,numelementsx) + repmat(((minix-ms:step:maxix-ms)-1)*size(image1_roi, 1), numelementsy,1))';
		s0B = permute(s0B(:), [2 3 1]);
		s1B = repmat((1:interrogationarea)',1,interrogationarea) + repmat(((1:interrogationarea)-1)*size(image1_roi, 1),interrogationarea,1);
		ss1B = repmat(s1B, [1, 1, size(s0B,3)])+repmat(s0B, [interrogationarea, interrogationarea, 1]);
		image1_cutB = image1_roi(ss1B);
		image2_cutB = image2_roi(ss1B);
		if do_pad==1 && passes == 1
			%subtract mean to avoid high frequencies at border of correlation:
			image1_cutB=image1_cutB-mean(image1_cutB,[1 2]);
			image2_cutB=image2_cutB-mean(image2_cutB,[1 2]);
			% padding (faster than padarray) to get the linear correlation:
			image1_cutB=[image1_cutB zeros(interrogationarea,interrogationarea-1,size(image1_cutB,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image1_cutB,3))];
			image2_cutB=[image2_cutB zeros(interrogationarea,interrogationarea-1,size(image2_cutB,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image2_cutB,3))];
		end
		result_convB = fftshift(fftshift(real(ifft2(conj(fft2(image1_cutB)).*fft2(image2_cutB))), 1), 2);
		if do_pad==1 && passes == 1
			%cropping of correlation matrix:
			result_convB =result_convB((interrogationarea/2):(3*interrogationarea/2)-1,(interrogationarea/2):(3*interrogationarea/2)-1,:);
		end

		%Shift right bot
		s0C = (repmat((miniy+ms:step:maxiy+ms)'-1, 1,numelementsx) + repmat(((minix+ms:step:maxix+ms)-1)*size(image1_roi, 1), numelementsy,1))';
		s0C = permute(s0C(:), [2 3 1]);
		s1C = repmat((1:interrogationarea)',1,interrogationarea) + repmat(((1:interrogationarea)-1)*size(image1_roi, 1),interrogationarea,1);
		ss1C = repmat(s1C, [1, 1, size(s0C,3)])+repmat(s0C, [interrogationarea, interrogationarea, 1]);
		image1_cutC = image1_roi(ss1C);
		image2_cutC = image2_roi(ss1C);
		if do_pad==1 && passes == 1
			%subtract mean to avoid high frequencies at border of correlation:
			image1_cutC=image1_cutC-mean(image1_cutC,[1 2]);
			image2_cutC=image2_cutC-mean(image2_cutC,[1 2]);
			% padding (faster than padarray) to get the linear correlation:
			image1_cutC=[image1_cutC zeros(interrogationarea,interrogationarea-1,size(image1_cutC,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image1_cutC,3))];
			image2_cutC=[image2_cutC zeros(interrogationarea,interrogationarea-1,size(image2_cutC,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image2_cutC,3))];
		end
		result_convC = fftshift(fftshift(real(ifft2(conj(fft2(image1_cutC)).*fft2(image2_cutC))), 1), 2);
		if do_pad==1 && passes == 1
			%cropping of correlation matrix:
			result_convC =result_convC((interrogationarea/2):(3*interrogationarea/2)-1,(interrogationarea/2):(3*interrogationarea/2)-1,:);
		end

		%Shift left top
		s0D = (repmat((miniy-ms:step:maxiy-ms)'-1, 1,numelementsx) + repmat(((minix-ms:step:maxix-ms)-1)*size(image1_roi, 1), numelementsy,1))';
		s0D = permute(s0D(:), [2 3 1]);
		s1D = repmat((1:interrogationarea)',1,interrogationarea) + repmat(((1:interrogationarea)-1)*size(image1_roi, 1),interrogationarea,1);
		ss1D = repmat(s1D, [1, 1, size(s0D,3)])+repmat(s0D, [interrogationarea, interrogationarea, 1]);
		image1_cutD = image1_roi(ss1D);
		image2_cutD = image2_roi(ss1D);

		if do_pad==1 && passes == 1
			%subtract mean to avoid high frequencies at border of correlation:
			image1_cutD=image1_cutD-mean(image1_cutD,[1 2]);
			image2_cutD=image2_cutD-mean(image2_cutD,[1 2]);
			% padding (faster than padarray) to get the linear correlation:
			image1_cutD=[image1_cutD zeros(interrogationarea,interrogationarea-1,size(image1_cutD,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image1_cutD,3))];
			image2_cutD=[image2_cutD zeros(interrogationarea,interrogationarea-1,size(image2_cutD,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image2_cutD,3))];
		end
		result_convD = fftshift(fftshift(real(ifft2(conj(fft2(image1_cutD)).*fft2(image2_cutD))), 1), 2);
		if do_pad==1 && passes == 1
			%cropping of correlation matrix:
			result_convD =result_convD((interrogationarea/2):(3*interrogationarea/2)-1,(interrogationarea/2):(3*interrogationarea/2)-1,:);
		end

		%Shift right top
		s0E = (repmat((miniy-ms:step:maxiy-ms)'-1, 1,numelementsx) + repmat(((minix+ms:step:maxix+ms)-1)*size(image1_roi, 1), numelementsy,1))';
		s0E = permute(s0E(:), [2 3 1]);
		s1E = repmat((1:interrogationarea)',1,interrogationarea) + repmat(((1:interrogationarea)-1)*size(image1_roi, 1),interrogationarea,1);
		ss1E = repmat(s1E, [1, 1, size(s0E,3)])+repmat(s0E, [interrogationarea, interrogationarea, 1]);
		image1_cutE = image1_roi(ss1E);
		image2_cutE = image2_roi(ss1E);
		if do_pad==1 && passes == 1
			%subtract mean to avoid high frequencies at border of correlation:
			image1_cutE=image1_cutE-mean(image1_cutE,[1 2]);
			image2_cutE=image2_cutE-mean(image2_cutE,[1 2]);
			% padding (faster than padarray) to get the linear correlation:
			image1_cutE=[image1_cutE zeros(interrogationarea,interrogationarea-1,size(image1_cutE,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image1_cutE,3))];
			image2_cutE=[image2_cutE zeros(interrogationarea,interrogationarea-1,size(image2_cutE,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image2_cutE,3))];
		end
		result_convE = fftshift(fftshift(real(ifft2(conj(fft2(image1_cutE)).*fft2(image2_cutE))), 1), 2);
		if do_pad==1 && passes == 1
			%cropping of correlation matrix:
			result_convE =result_convE((interrogationarea/2):(3*interrogationarea/2)-1,(interrogationarea/2):(3*interrogationarea/2)-1,:);
		end
		result_conv=result_conv.*result_convB.*result_convC.*result_convD.*result_convE;
	end

	if mask_auto == 1
		%das zentrum der Matrize (3x3) mit dem mittelwert ersetzen = Keine Autokorrelation
		%MARKER
		h = fspecial('gaussian', 3, 1.5);
		h=h/h(2,2);
		h=1-h;
		%h=repmat(h,1,1,size(result_conv,3));
		h=repmat(h,[1,1,size(result_conv,3)]);
		h=h.*result_conv((interrogationarea/2)+SubPixOffset-1:(interrogationarea/2)+SubPixOffset+1,(interrogationarea/2)+SubPixOffset-1:(interrogationarea/2)+SubPixOffset+1,:);
		result_conv((interrogationarea/2)+SubPixOffset-1:(interrogationarea/2)+SubPixOffset+1,(interrogationarea/2)+SubPixOffset-1:(interrogationarea/2)+SubPixOffset+1,:)=h;
	end

    %% Apply mask — piv_FFTensemble lines 287-288
    ii=find(mask(ss1(round(interrogationarea/2+1),round(interrogationarea/2+1),:)));
    result_conv(:,:,ii)=0;

    %% Accumulate — piv_FFTensemble lines 290-295
    try
        result_conv_ensemble=result_conv_ensemble+result_conv;
    catch
        result_conv_ensemble=zeros(size(result_conv));
        result_conv_ensemble=result_conv_ensemble+result_conv;
    end

    %% Correlation map (pass 1 only) — piv_FFTensemble lines 348-357
    if passes==1
        for cor_i=1:size(image1_cut,3)
            correlation_map(cor_i)=correlation_map(cor_i)+corr2(image1_cut(:,:,cor_i),image2_cut(:,:,cor_i));
        end
        corr_map_cnt=corr_map_cnt+1;
    end

end % frame loop pass 1
fprintf('Done.\n');

%% Pass-1 peakfinding — piv_FFTensemble line 362
if use_gpu
    result_conv_ensemble = gather(result_conv_ensemble);
    mask = gather(mask);
    ss1 = gather(ss1);
end
[xtable,ytable,utable,vtable]=peakfinding(result_conv_ensemble,mask,interrogationarea,minix,step,maxix,miniy,maxiy,SubPixOffset,ss1,subpixfinder);

%% Multipass refinement — piv_FFTensemble lines 363-806
for multipass=1:passes-1
    if multipass==1; interrogationarea=round(int2/2)*2; end % piv_FFTensemble lines 365-373
    if multipass==2; interrogationarea=round(int3/2)*2; end
    if multipass==3; interrogationarea=round(int4/2)*2; end
    % [DEVIATION] piv_FFTensemble: step=interrogationarea/2 (always 50%).
    % Here: step from window_sizes, allowing custom overlap.
    step=ss(multipass+1);

    result_conv_ensemble=zeros(interrogationarea,interrogationarea); % piv_FFTensemble line 374

    %% Predictor validation/smoothing — piv_FFTensemble lines 468-512
    % [DEVIATION] piv_FFTensemble puts this inside the per-frame loop. Since
    % utable/vtable don't change within a pass, the result is the same each
    % iteration there. We move it outside so it runs ONCE per pass, preventing
    % progressive over-smoothing across n_pairs iterations.
    if use_gpu
        utable = gather(utable); vtable = gather(vtable);
    end
    utable_orig=utable; vtable_orig=vtable;
    [utable,vtable]=postproc.PIVlab_postproc(utable,vtable,[],[],[],1,4,1,1.5);
    utable=misc.inpaint_nans(utable,4); vtable=misc.inpaint_nans(vtable,4);
    if multipass<passes-1
        utable=misc.smoothn(utable,0.6); vtable=misc.smoothn(vtable,0.6);
    else
        utable=misc.smoothn(utable); vtable=misc.smoothn(vtable);
    end

    %% Restore gen_image, recompute grid
    image1_size_1 = heightroi+1;
    image1_size_2 = widthroi+1;

    miniy=1+(ceil(interrogationarea/2));
    minix=1+(ceil(interrogationarea/2));
    maxiy=step*(floor(image1_size_1/step))-(interrogationarea-1)+(ceil(interrogationarea/2));
    maxix=step*(floor(image1_size_2/step))-(interrogationarea-1)+(ceil(interrogationarea/2));

    numelementsy=floor((maxiy-miniy)/step+1);
    numelementsx=floor((maxix-minix)/step+1);

    LAy=miniy; LAx=minix;
    LUy=image1_size_1-maxiy; LUx=image1_size_2-maxix;
    shift4centery=round((LUy-LAy)/2); shift4centerx=round((LUx-LAx)/2);
    if shift4centery<0; shift4centery=0; end
    if shift4centerx<0; shift4centerx=0; end
    miniy=miniy+shift4centery; minix=minix+shift4centerx;
    maxix=maxix+shift4centerx; maxiy=maxiy+shift4centery;

    if (rem(interrogationarea,2)==0); SubPixOffset=1; else; SubPixOffset=0.5; end

    xtable_old=xtable; ytable_old=ytable;
    typevector=ones(numelementsy,numelementsx);
    xtable=repmat((minix:step:maxix),numelementsy,1)+interrogationarea/2;
    ytable=repmat((miniy:step:maxiy)',1,numelementsx)+interrogationarea/2;

    %% Predictor interpolation MUST happen on CPU with *spline to handle edge extrapolation safely
    if use_gpu
        utable = gather(utable); vtable = gather(vtable);
        xtable_old = gather(xtable_old); ytable_old = gather(ytable_old);
    end
    utable=interp2(xtable_old,ytable_old,utable,xtable,ytable,'*spline');
    vtable=interp2(xtable_old,ytable_old,vtable,xtable,ytable,'*spline');
    
    %% Pixel-resolution predictor field
    utable_1=padarray(utable,[1,1],'replicate');
    vtable_1=padarray(vtable,[1,1],'replicate');
    firstlinex=xtable(1,:);
    firstlinex_intp=interp1(1:size(firstlinex,2),firstlinex,0:size(firstlinex,2)+1,'linear','extrap');
    xtable_1=repmat(firstlinex_intp,size(xtable,1)+2,1);
    firstliney=ytable(:,1);
    firstliney_intp=interp1(1:size(firstliney,1),firstliney,0:size(firstliney,1)+1,'linear','extrap')';
    ytable_1=repmat(firstliney_intp,1,size(ytable,2)+2);
    X=xtable_1; Y=ytable_1; U=utable_1; V=vtable_1;
    X1=X(1,1):1:X(1,end)-1;
    Y1=(Y(1,1):1:Y(end,1)-1)';
    X1=repmat(X1,size(Y1,1),1);
    Y1=repmat(Y1,1,size(X1,2));
    U1=interp2(X,Y,U,X1,Y1,'*linear');
    V1=interp2(X,Y,V,X1,Y1,'*linear');

    image2_x_grid = 1:(image1_size_2 + 2*ceil(interrogationarea/2));
    image2_y_grid = (1:(image1_size_1 + 2*ceil(interrogationarea/2)))';

    if use_gpu
        xtable = gpuArray(xtable); ytable = gpuArray(ytable);
        U1 = gpuArray(U1); V1 = gpuArray(V1);
        X1 = gpuArray(X1); Y1 = gpuArray(Y1);
        xtable_1 = gpuArray(xtable_1); ytable_1 = gpuArray(ytable_1);
        image2_x_grid = gpuArray(image2_x_grid); image2_y_grid = gpuArray(image2_y_grid);
    end

    if use_gpu
        fprintf('Pass %d/%d: Processing %d pairs on GPU... ', multipass+1, passes, n_pairs);
    else
        fprintf('Pass %d/%d: Processing %d pairs on CPU... ', multipass+1, passes, n_pairs);
    end

    for k=1:n_pairs % [DEVIATION] consecutive pairs
        if mod(k, 10) == 0; fprintf('.'); if use_gpu; wait(gpuDevice); end; end

        %% [DEVIATION] Load images from imgstack
        image1_roi=double(imgstack(yroi:yroi+heightroi,xroi:xroi+widthroi,k));
        image2_roi=double(imgstack(yroi:yroi+heightroi,xroi:xroi+widthroi,k+1));
        mask_inpt_roi=mask_inpt(yroi:yroi+heightroi,xroi:xroi+widthroi);
        mask=mask_inpt_roi;

        %% Pad
        image1_roi=padarray(image1_roi,[ceil(interrogationarea/2) ceil(interrogationarea/2)],min(min(image1_roi)));
        image2_roi=padarray(image2_roi,[ceil(interrogationarea/2) ceil(interrogationarea/2)],min(min(image1_roi)));
        mask      =padarray(mask,      [ceil(interrogationarea/2) ceil(interrogationarea/2)],0);

        %% Image deformation
        image2_crop_i1=interp2(image2_x_grid,image2_y_grid,double(image2_roi),X1+U1,Y1+V1,imdeform,min(image2_roi(:)));

        xb=find(X1(1,:)==xtable_1(1,1)); yb=find(Y1(:,1)==ytable_1(1,1));

        %% ss1 (image1) and ss2 (image2_crop) index arrays — piv_FFTensemble lines 609-620
        s0=(repmat((miniy:step:maxiy)'-1,1,numelementsx)+repmat(((minix:step:maxix)-1)*size(image1_roi,1),numelementsy,1))';
        s0=permute(s0(:),[2 3 1]);
        s1=repmat((1:interrogationarea)',1,interrogationarea)+repmat(((1:interrogationarea)-1)*size(image1_roi,1),interrogationarea,1);
        ss1=repmat(s1,[1,1,size(s0,3)])+repmat(s0,[interrogationarea,interrogationarea,1]);

        s0b=(repmat(yb-step+step*(1:numelementsy)'-1,1,numelementsx)+repmat((xb-step+step*(1:numelementsx)-1)*size(image2_crop_i1,1),numelementsy,1))';
        s0b=permute(s0b(:),[2 3 1])-s0b(1);
        s2=repmat((1:2*step)',1,2*step)+repmat(((1:2*step)-1)*size(image2_crop_i1,1),2*step,1);
        ss2=repmat(s2,[1,1,size(s0b,3)])+repmat(s0b,[interrogationarea,interrogationarea,1]);

        image1_cut=image1_roi(ss1);
        image2_cut=image2_crop_i1(ss2);

        %% do_pad on last pass — piv_FFTensemble lines 621-635
        is_last=(multipass==passes-1);
        if do_pad==1 && is_last
            image1_cut=image1_cut-mean(image1_cut,[1 2]);
            image2_cut=image2_cut-mean(image2_cut,[1 2]);
            image1_cut=[image1_cut zeros(interrogationarea,interrogationarea-1,size(image1_cut,3));zeros(interrogationarea-1,2*interrogationarea-1,size(image1_cut,3))];
            image2_cut=[image2_cut zeros(interrogationarea,interrogationarea-1,size(image2_cut,3));zeros(interrogationarea-1,2*interrogationarea-1,size(image2_cut,3))];
        end
        result_conv=fftshift(fftshift(real(ifft2(conj(fft2(image1_cut)).*fft2(image2_cut))),1),2);
        if do_pad==1 && is_last
            result_conv=result_conv((interrogationarea/2):(3*interrogationarea/2)-1,(interrogationarea/2):(3*interrogationarea/2)-1,:);
        end

			%% repeated correlation
			if repeat == 1 && multipass==passes-1
				ms=round(step/4); %multishift parameter so groß wie viertel int window

				%Shift left bot
				image2_crop_i1 = interp2(image2_x_grid,image2_y_grid,double(image2_roi),X1+U1-ms,Y1+V1+ms,imdeform,min(image2_roi(:))); %linear is 3x faster and looks ok...
				xb = find(X1(1,:) == xtable_1(1,1));
				yb = find(Y1(:,1) == ytable_1(1,1));
				s0 = (repmat((miniy+ms:step:maxiy+ms)'-1, 1,numelementsx) + repmat(((minix-ms:step:maxix-ms)-1)*size(image1_roi, 1), numelementsy,1))';
				s0 = permute(s0(:), [2 3 1]);
				s1 = repmat((1:interrogationarea)',1,interrogationarea) + repmat(((1:interrogationarea)-1)*size(image1_roi, 1),interrogationarea,1);
				ss1 = repmat(s1, [1, 1, size(s0,3)]) + repmat(s0, [interrogationarea, interrogationarea, 1]);
				s0 = (repmat(yb-step+step*(1:numelementsy)'-1, 1,numelementsx) + repmat((xb-step+step*(1:numelementsx)-1)*size(image2_crop_i1, 1), numelementsy,1))';
				s0 = permute(s0(:), [2 3 1]) - s0(1);
				s2 = repmat((1:2*step)',1,2*step) + repmat(((1:2*step)-1)*size(image2_crop_i1, 1),2*step,1);
				ss2 = repmat(s2, [1, 1, size(s0,3)]) + repmat(s0, [interrogationarea, interrogationarea, 1]);
				image1_cut = image1_roi(ss1);
				image2_cut = image2_crop_i1(ss2);
				if do_pad==1 && multipass==passes-1
					%subtract mean to avoid high frequencies at border of correlation:

					image1_cut=image1_cut-mean(image1_cut,[1 2]);
					image2_cut=image2_cut-mean(image2_cut,[1 2]);

					% padding (faster than padarray) to get the linear correlation:
					image1_cut=[image1_cut zeros(interrogationarea,interrogationarea-1,size(image1_cut,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image1_cut,3))];
					image2_cut=[image2_cut zeros(interrogationarea,interrogationarea-1,size(image2_cut,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image2_cut,3))];
				end
				result_convB = fftshift(fftshift(real(ifft2(conj(fft2(image1_cut)).*fft2(image2_cut))), 1), 2);
				if do_pad==1 && multipass==passes-1
					%cropping of correlation matrix:
					result_convB =result_convB((interrogationarea/2):(3*interrogationarea/2)-1,(interrogationarea/2):(3*interrogationarea/2)-1,:);
				end
				%Shift right bot
				image2_crop_i1 = interp2(image2_x_grid,image2_y_grid,double(image2_roi),X1+U1+ms,Y1+V1+ms,imdeform,min(image2_roi(:))); %linear is 3x faster and looks ok...
				xb = find(X1(1,:) == xtable_1(1,1));
				yb = find(Y1(:,1) == ytable_1(1,1));
				s0 = (repmat((miniy+ms:step:maxiy+ms)'-1, 1,numelementsx) + repmat(((minix+ms:step:maxix+ms)-1)*size(image1_roi, 1), numelementsy,1))';
				s0 = permute(s0(:), [2 3 1]);
				s1 = repmat((1:interrogationarea)',1,interrogationarea) + repmat(((1:interrogationarea)-1)*size(image1_roi, 1),interrogationarea,1);
				ss1 = repmat(s1, [1, 1, size(s0,3)]) + repmat(s0, [interrogationarea, interrogationarea, 1]);
				s0 = (repmat(yb-step+step*(1:numelementsy)'-1, 1,numelementsx) + repmat((xb-step+step*(1:numelementsx)-1)*size(image2_crop_i1, 1), numelementsy,1))';
				s0 = permute(s0(:), [2 3 1]) - s0(1);
				s2 = repmat((1:2*step)',1,2*step) + repmat(((1:2*step)-1)*size(image2_crop_i1, 1),2*step,1);
				ss2 = repmat(s2, [1, 1, size(s0,3)]) + repmat(s0, [interrogationarea, interrogationarea, 1]);
				image1_cut = image1_roi(ss1);
				image2_cut = image2_crop_i1(ss2);
				if do_pad==1 && multipass==passes-1
					%subtract mean to avoid high frequencies at border of correlation:

					image1_cut=image1_cut-mean(image1_cut,[1 2]);
					image2_cut=image2_cut-mean(image2_cut,[1 2]);

					% padding (faster than padarray) to get the linear correlation:
					image1_cut=[image1_cut zeros(interrogationarea,interrogationarea-1,size(image1_cut,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image1_cut,3))];
					image2_cut=[image2_cut zeros(interrogationarea,interrogationarea-1,size(image2_cut,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image2_cut,3))];
				end
				result_convC = fftshift(fftshift(real(ifft2(conj(fft2(image1_cut)).*fft2(image2_cut))), 1), 2);
				if do_pad==1 && multipass==passes-1
					%cropping of correlation matrix:
					result_convC =result_convC((interrogationarea/2):(3*interrogationarea/2)-1,(interrogationarea/2):(3*interrogationarea/2)-1,:);
				end
				%Shift left top
				image2_crop_i1 = interp2(image2_x_grid,image2_y_grid,double(image2_roi),X1+U1-ms,Y1+V1-ms,imdeform,min(image2_roi(:))); %linear is 3x faster and looks ok...
				xb = find(X1(1,:) == xtable_1(1,1));
				yb = find(Y1(:,1) == ytable_1(1,1));
				s0 = (repmat((miniy-ms:step:maxiy-ms)'-1, 1,numelementsx) + repmat(((minix-ms:step:maxix-ms)-1)*size(image1_roi, 1), numelementsy,1))';
				s0 = permute(s0(:), [2 3 1]);
				s1 = repmat((1:interrogationarea)',1,interrogationarea) + repmat(((1:interrogationarea)-1)*size(image1_roi, 1),interrogationarea,1);
				ss1 = repmat(s1, [1, 1, size(s0,3)]) + repmat(s0, [interrogationarea, interrogationarea, 1]);
				s0 = (repmat(yb-step+step*(1:numelementsy)'-1, 1,numelementsx) + repmat((xb-step+step*(1:numelementsx)-1)*size(image2_crop_i1, 1), numelementsy,1))';
				s0 = permute(s0(:), [2 3 1]) - s0(1);
				s2 = repmat((1:2*step)',1,2*step) + repmat(((1:2*step)-1)*size(image2_crop_i1, 1),2*step,1);
				ss2 = repmat(s2, [1, 1, size(s0,3)]) + repmat(s0, [interrogationarea, interrogationarea, 1]);
				image1_cut = image1_roi(ss1);
				image2_cut = image2_crop_i1(ss2);
				if do_pad==1 && multipass==passes-1
					%subtract mean to avoid high frequencies at border of correlation:

					image1_cut=image1_cut-mean(image1_cut,[1 2]);
					image2_cut=image2_cut-mean(image2_cut,[1 2]);

					% padding (faster than padarray) to get the linear correlation:
					image1_cut=[image1_cut zeros(interrogationarea,interrogationarea-1,size(image1_cut,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image1_cut,3))];
					image2_cut=[image2_cut zeros(interrogationarea,interrogationarea-1,size(image2_cut,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image2_cut,3))];
				end
				result_convD = fftshift(fftshift(real(ifft2(conj(fft2(image1_cut)).*fft2(image2_cut))), 1), 2);
				if do_pad==1 && multipass==passes-1
					%cropping of correlation matrix:
					result_convD =result_convD((interrogationarea/2):(3*interrogationarea/2)-1,(interrogationarea/2):(3*interrogationarea/2)-1,:);
				end
				%Shift right top
				image2_crop_i1 = interp2(image2_x_grid,image2_y_grid,double(image2_roi),X1+U1+ms,Y1+V1-ms,imdeform,min(image2_roi(:))); %linear is 3x faster and looks ok...
				xb = find(X1(1,:) == xtable_1(1,1));
				yb = find(Y1(:,1) == ytable_1(1,1));
				s0 = (repmat((miniy-ms:step:maxiy-ms)'-1, 1,numelementsx) + repmat(((minix+ms:step:maxix+ms)-1)*size(image1_roi, 1), numelementsy,1))';
				s0 = permute(s0(:), [2 3 1]);
				s1 = repmat((1:interrogationarea)',1,interrogationarea) + repmat(((1:interrogationarea)-1)*size(image1_roi, 1),interrogationarea,1);
				ss1 = repmat(s1, [1, 1, size(s0,3)]) + repmat(s0, [interrogationarea, interrogationarea, 1]);
				s0 = (repmat(yb-step+step*(1:numelementsy)'-1, 1,numelementsx) + repmat((xb-step+step*(1:numelementsx)-1)*size(image2_crop_i1, 1), numelementsy,1))';
				s0 = permute(s0(:), [2 3 1]) - s0(1);
				s2 = repmat((1:2*step)',1,2*step) + repmat(((1:2*step)-1)*size(image2_crop_i1, 1),2*step,1);
				ss2 = repmat(s2, [1, 1, size(s0,3)]) + repmat(s0, [interrogationarea, interrogationarea, 1]);
				image1_cut = image1_roi(ss1);
				image2_cut = image2_crop_i1(ss2);
				if do_pad==1 && multipass==passes-1
					%subtract mean to avoid high frequencies at border of correlation:

					image1_cut=image1_cut-mean(image1_cut,[1 2]);
					image2_cut=image2_cut-mean(image2_cut,[1 2]);

					% padding (faster than padarray) to get the linear correlation:
					image1_cut=[image1_cut zeros(interrogationarea,interrogationarea-1,size(image1_cut,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image1_cut,3))];
					image2_cut=[image2_cut zeros(interrogationarea,interrogationarea-1,size(image2_cut,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image2_cut,3))];
				end
				result_convE = fftshift(fftshift(real(ifft2(conj(fft2(image1_cut)).*fft2(image2_cut))), 1), 2);
				if do_pad==1 && multipass==passes-1
					%cropping of correlation matrix:
					result_convE =result_convE((interrogationarea/2):(3*interrogationarea/2)-1,(interrogationarea/2):(3*interrogationarea/2)-1,:);
				end
				result_conv=result_conv.*result_convB.*result_convC.*result_convD.*result_convE;
			end


			if mask_auto == 1
				%limit peak search arena....
				emptymatrix=zeros(size(result_conv,1),size(result_conv,2),size(result_conv,3));
				%emptymatrix=emptymatrix+0.1;
				if interrogationarea > 8 % masking central peak will not work for extrmely small interrogation areas. And it also doesn't make sense.
					sizeones=4;
					%h = fspecial('gaussian', sizeones*2+1,1);
					h=fspecial('disk',4);
					h=h/max(max(h));
					%h=repmat(h,1,1,size(result_conv,3));
					h=repmat(h,[1,1,size(result_conv,3)]);
					emptymatrix((interrogationarea/2)+SubPixOffset-sizeones:(interrogationarea/2)+SubPixOffset+sizeones,(interrogationarea/2)+SubPixOffset-sizeones:(interrogationarea/2)+SubPixOffset+sizeones,:)=h;
					result_conv = result_conv .* emptymatrix;
				else
					disp('All interrogation areas must be larger than 8 pixels for disabling auto correlation successfully.')
				end
			end

        %% Apply mask — piv_FFTensemble lines 780-781
        ii=find(mask(ss1(round(interrogationarea/2+1),round(interrogationarea/2+1),:)));
        result_conv(:,:,ii)=0;

        %% Accumulate — piv_FFTensemble lines 783-788
        try
            result_conv_ensemble=result_conv_ensemble+result_conv;
        catch
            result_conv_ensemble=zeros(size(result_conv));
            result_conv_ensemble=result_conv_ensemble+result_conv;
        end

        %% Correlation map on last pass — piv_FFTensemble lines 790-801
        if is_last
            if k==1; correlation_map=zeros(size(typevector)); corr_map_cnt=0; end
            for cor_i=1:size(image1_cut,3)
                correlation_map(cor_i)=correlation_map(cor_i)+corr2(image1_cut(:,:,cor_i),image2_cut(:,:,cor_i));
            end
            corr_map_cnt=corr_map_cnt+1;
        end

    end % frame loop multipass
fprintf('Done.\n');

    %% Peakfinding and accumulate displacement — piv_FFTensemble lines 803-805
    if use_gpu
        result_conv_ensemble = gather(result_conv_ensemble);
        ss1 = gather(ss1);
    end
    [xtable,ytable,utable2,vtable2]=peakfinding(result_conv_ensemble,[],interrogationarea,minix,step,maxix,miniy,maxiy,SubPixOffset,ss1,subpixfinder);
    utable=utable+utable2;
    vtable=vtable+vtable2;

end % multipass loop

%% Apply average mask — piv_FFTensemble lines 807-828
nrx=0; nrxreal=0; nry=0;
average_mask_roi_pad=padarray(average_mask_roi,[ceil(interrogationarea/2) ceil(interrogationarea/2)],0);
for jmask=miniy:step:maxiy
    nry=nry+1;
    for imask=minix:step:maxix
        nrx=nrx+1;
        if nrxreal<numelementsx; nrxreal=nrxreal+1; else; nrxreal=1; end
        if average_mask_roi_pad(round(jmask+interrogationarea/2),round(imask+interrogationarea/2))>=n_pairs
            typevector(nry,nrxreal)=0;
        end
    end
end

%% Coordinate correction to global frame — piv_FFTensemble lines 829-833
xtable=xtable-ceil(interrogationarea/2);
ytable=ytable-ceil(interrogationarea/2);
xtable=xtable+xroi;
ytable=ytable+yroi;

%% Normalise correlation map — piv_FFTensemble line 852
correlation_map=permute(reshape(correlation_map,[size(xtable')]),[2 1 3])/corr_map_cnt;
correlation_map(typevector==0)=0;

end % piv_ensemble

%%-------------------------------------------------------------------------
function [xtable,ytable,utable,vtable] = peakfinding(result_conv_ensemble, mask, interrogationarea, minix, step, maxix, miniy, maxiy, SubPixOffset, ss1, subpixfinder)
minres = permute(repmat(squeeze(min(min(result_conv_ensemble))), [1, size(result_conv_ensemble,1), size(result_conv_ensemble,2)]), [2 3 1]);
deltares = permute(repmat(squeeze(max(max(result_conv_ensemble))-min(min(result_conv_ensemble))), [1, size(result_conv_ensemble,1), size(result_conv_ensemble,2)]), [2 3 1]);
result_conv_ensemble = ((result_conv_ensemble-minres)./deltares)*255;

%apply mask
if isempty(mask)==0
	ii = find(mask(ss1(round(interrogationarea/2+1), round(interrogationarea/2+1), :)));
	result_conv_ensemble(:,:,ii) = 0;
end

[y, x, z] = ind2sub(size(result_conv_ensemble), find(result_conv_ensemble==255));

% we need only one peak from each couple pictures
[z1, zi] = sort(z);
dz1 = [z1(1); diff(z1)];
i0 = find(dz1~=0);
x1 = x(zi(i0));
y1 = y(zi(i0));
z1 = z(zi(i0));

xtable = repmat((minix:step:maxix)+interrogationarea/2, length(miniy:step:maxiy), 1);
ytable = repmat(((miniy:step:maxiy)+interrogationarea/2)', 1, length(minix:step:maxix));

if subpixfinder==1
	[vector] = SUBPIXGAUSS(result_conv_ensemble, interrogationarea, x1, y1, z1, SubPixOffset);
elseif subpixfinder==2
	[vector] = SUBPIX2DGAUSS(result_conv_ensemble, interrogationarea, x1, y1, z1, SubPixOffset);
end
vector = permute(reshape(vector, [size(xtable') 2]), [2 1 3]);
utable = vector(:,:,1);
vtable = vector(:,:,2);
end

function vector = SUBPIXGAUSS(result_conv, interrogationarea, x, y, z, SubPixOffset)
    xi = find(~((x <= size(result_conv,2)-1) & (y <= size(result_conv,1)-1) & (x >= 2) & (y >= 2)));
    x(xi) = []; y(xi) = []; z(xi) = [];
    xmax = size(result_conv, 2);
    vector = NaN(size(result_conv,3), 2);
    if numel(x) ~= 0
        ip = sub2ind(size(result_conv), y, x, z);
        f0 = log(result_conv(ip));
        f1 = log(result_conv(ip-1));   f2 = log(result_conv(ip+1));
        peaky = y + (f1-f2) ./ (2*f1 - 4*f0 + 2*f2);
        f1 = log(result_conv(ip-xmax)); f2 = log(result_conv(ip+xmax));
        peakx = x + (f1-f2) ./ (2*f1 - 4*f0 + 2*f2);
        vector(z,:) = [peakx - interrogationarea/2 - SubPixOffset, ...
                       peaky - interrogationarea/2 - SubPixOffset];
    end
end


function vector = SUBPIX2DGAUSS(result_conv, interrogationarea, x, y, z, SubPixOffset)
    xi = find(~((x <= size(result_conv,2)-1) & (y <= size(result_conv,1)-1) & (x >= 2) & (y >= 2)));
    x(xi) = []; y(xi) = []; z(xi) = [];
    xmax = size(result_conv, 2);
    vector = NaN(size(result_conv,3), 2);
        if numel(x) ~= 0
            c10 = zeros(3,3,length(z)); c01=c10; c11=c10; c20=c10; c02=c10;
            ip = sub2ind(size(result_conv), y, x, z);
            for i = -1:1
                for j = -1:1
                    c10(j+2,i+2,:) = i*log(result_conv(ip+xmax*i+j));
                    c01(j+2,i+2,:) = j*log(result_conv(ip+xmax*i+j));
                    c11(j+2,i+2,:) = i*j*log(result_conv(ip+xmax*i+j));
                    c20(j+2,i+2,:) = (3*i^2-2)*log(result_conv(ip+xmax*i+j));
                    c02(j+2,i+2,:) = (3*j^2-2)*log(result_conv(ip+xmax*i+j));
                end
            end
            c10=(1/6)*sum(sum(c10)); c01=(1/6)*sum(sum(c01));
            c11=(1/4)*sum(sum(c11)); c20=(1/6)*sum(sum(c20)); c02=(1/6)*sum(sum(c02));
            deltax = squeeze((c11.*c01 - 2*c10.*c02) ./ (4*c20.*c02 - c11.^2));
            deltay = squeeze((c11.*c10 - 2*c01.*c20) ./ (4*c20.*c02 - c11.^2));
            peakx = x + deltax; peaky = y + deltay;
            vector(z,:) = [peakx - interrogationarea/2 - SubPixOffset, ...
                           peaky - interrogationarea/2 - SubPixOffset];
        end
end










