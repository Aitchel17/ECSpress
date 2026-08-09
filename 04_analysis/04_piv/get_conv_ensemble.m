function get_conv_ensemble(interrogationarea)

    %% pre-processing is done in this function
result_conv_ensemble = zeros(interrogationarea,interrogationarea); % prepare empty result_conv

num_slice = numel(video_frame_selection);
total_analyses_amount = num_slice / 2 * passes;
from_total = 0;
tic

    for ensemble_i1=1:2:num_slice

	image1 = read(filepath,video_frame_selection(ensemble_i1));
	image2 = read(filepath,video_frame_selection(ensemble_i1+1));

	if size(image1,3)>1
		image1=uint8(mean(image1,3));
		image2=uint8(mean(image2,3));
		%disp('Warning: To optimize speed, your images should be grayscale, 8 bit!')
	end

	%if autolimit == 1 %if autolimit is desired: do autolimit for each image seperately
	if size(image1,3)>1
		stretcher = stretchlim(rgb2gray(image1));
	else
		stretcher = stretchlim(image1);
	end
	minintens1 = stretcher(1);
	maxintens1 = stretcher(2);
	if size(image2,3)>1
		stretcher = stretchlim(rgb2gray(image2));
	else
		stretcher = stretchlim(image2);
	end
	minintens2 = stretcher(1);
	maxintens2 = stretcher(2);

	%% calculate the average mask that will be applied in the very end. It will remove all vectors where 100% of the input images have been masked.
	%prepare a matrix for calculating the average mask of all images
	if ensemble_i1==1
		average_mask=zeros(size(converted_mask{1,1}));
	end
	mask_inpt = converted_mask{floor((ensemble_i1+1)/2),1};
	if numel(roi_inpt)>0
		xroi = roi_inpt(1);
		yroi = roi_inpt(2);
		widthroi = roi_inpt(3);
		heightroi = roi_inpt(4);
		image1_roi = double(image1(yroi:yroi+heightroi,xroi:xroi+widthroi));
		image2_roi = double(image2(yroi:yroi+heightroi,xroi:xroi+widthroi));
		mask_inpt_roi = mask_inpt(yroi:yroi+heightroi,xroi:xroi+widthroi);
		average_mask_roi = average_mask(yroi:yroi+heightroi,xroi:xroi+widthroi);
	else
		xroi=0;
		yroi=0;
		image1_roi=double(image1);
		image2_roi=double(image2);
		mask_inpt_roi = mask_inpt;
		average_mask_roi = average_mask;
	end
	mask=mask_inpt_roi;
	gen_image1_roi = image1_roi;
	gen_image2_roi = image2_roi;


	miniy=1+(ceil(interrogationarea/2));
	minix=1+(ceil(interrogationarea/2));
	maxiy=step*(floor(size(image1_roi,1)/step))-(interrogationarea-1)+(ceil(interrogationarea/2)); %statt size deltax von ROI nehmen
	maxix=step*(floor(size(image1_roi,2)/step))-(interrogationarea-1)+(ceil(interrogationarea/2));

	numelementsy=floor((maxiy-miniy)/step+1);
	numelementsx=floor((maxix-minix)/step+1);

	LAy=miniy;
	LAx=minix;
	LUy=size(image1_roi,1)-maxiy;
	LUx=size(image1_roi,2)-maxix;
	shift4centery=round((LUy-LAy)/2);
	shift4centerx=round((LUx-LAx)/2);
	if shift4centery<0 %shift4center will be negative if in the unshifted case the left border is bigger than the right border. the vectormatrix is hence not centered on the image. the matrix cannot be shifted more towards the left border because then image2_crop would have a negative index. The only way to center the matrix would be to remove a column of vectors on the right side. but then we weould have less data....
		shift4centery=0;
	end
	if shift4centerx<0 %shift4center will be negative if in the unshifted case the left border is bigger than the right border. the vectormatrix is hence not centered on the image. the matrix cannot be shifted more towards the left border because then image2_crop would have a negative index. The only way to center the matrix would be to remove a column of vectors on the right side. but then we weould have less data....
		shift4centerx=0;
	end
	miniy=miniy+shift4centery;
	minix=minix+shift4centerx;
	maxix=maxix+shift4centerx;
	maxiy=maxiy+shift4centery;

	image1_roi=padarray(image1_roi,[ceil(interrogationarea/2) ceil(interrogationarea/2)], min(min(image1_roi)));
	image2_roi=padarray(image2_roi,[ceil(interrogationarea/2) ceil(interrogationarea/2)], min(min(image1_roi)));
	mask=padarray(mask,[ceil(interrogationarea/2) ceil(interrogationarea/2)],0);

	if (rem(interrogationarea,2) == 0) %for the subpixel displacement measurement
		SubPixOffset=1;
	else
		SubPixOffset=0.5;
	end
	xtable=zeros(numelementsy,numelementsx);
	ytable=xtable; %#ok<*NASGU>
	utable=xtable;
	vtable=xtable;
	typevector=ones(numelementsy,numelementsx);

	%% MAINLOOP

	% divide images by small pictures
	% new index for image1_roi and image2_roi
	s0 = (repmat((miniy:step:maxiy)'-1, 1,numelementsx) + repmat(((minix:step:maxix)-1)*size(image1_roi, 1), numelementsy,1))';
	s0 = permute(s0(:), [2 3 1]);
	s1 = repmat((1:interrogationarea)',1,interrogationarea) + repmat(((1:interrogationarea)-1)*size(image1_roi, 1),interrogationarea,1);
	ss1 = repmat(s1, [1, 1, size(s0,3)])+repmat(s0, [interrogationarea, interrogationarea, 1]);

	image1_cut = image1_roi(ss1);
	image2_cut = image2_roi(ss1);

	if do_pad==1 && passes == 1 %only on first pass
		%subtract mean to avoid high frequencies at border of correlation:
		image1_cut=image1_cut-mean(image1_cut,[1 2]);
		image2_cut=image2_cut-mean(image2_cut,[1 2]);
		% padding (faster than padarray) to get the linear correlation:
		image1_cut=[image1_cut zeros(interrogationarea,interrogationarea-1,size(image1_cut,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image1_cut,3))];
		image2_cut=[image2_cut zeros(interrogationarea,interrogationarea-1,size(image2_cut,3)); zeros(interrogationarea-1,2*interrogationarea-1,size(image2_cut,3))];
	end
	%do fft2:

	result_conv = fftshift(fftshift(real(ifft2(conj(fft2(image1_cut)).*fft2(image2_cut))), 1), 2);
	if do_pad==1 && passes == 1
		%cropping of correlation matrix:
		result_conv =result_conv((interrogationarea/2):(3*interrogationarea/2)-1,(interrogationarea/2):(3*interrogationarea/2)-1,:);
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
	%apply mask
	ii = find(mask(ss1(round(interrogationarea/2+1), round(interrogationarea/2+1), :)));
	result_conv(:,:, ii) = 0;
	%average the correlation matrices
	try
		result_conv_ensemble=result_conv_ensemble+result_conv;
	catch % older matlab releases
		result_conv_ensemble = zeros(size(result_conv));
		result_conv_ensemble=result_conv_ensemble+result_conv;
	end

    fprintf('.');
	
	if passes==1 % only 1 pass selected, so correlation coefficient will be calculated in this (first & final) pass.
		if ensemble_i1==1 %first image pair
			correlation_map=zeros(size(typevector));
			corr_map_cnt=0;
		end
		for cor_i=1:size(image1_cut,3)
			correlation_map(cor_i)=correlation_map(cor_i) + corr2(image1_cut(:,:,cor_i),image2_cut(:,:,cor_i));
		end
		corr_map_cnt=corr_map_cnt+1;
	end
    end

end