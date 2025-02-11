clear ,clc, tic;
analysis = analyze('E:\OneDrive - The Pennsylvania State University\2023_ECSpress\01_primary_analysis\matlab_analysis');
%%





%%
analysis2 = analyze('E:\OneDrive - The Pennsylvania State University\2023_ECSpress\01_primary_analysis\matlab_analysis');

%% Extra parenchymal space

tmp.save_resolution = [str2double(analysis2.info.objpix(1:end-2)),str2double(analysis2.info.objpix(1:end-2)),1/str2double(analysis2.info.savefps)]; % [x,y,z resolution um, sec]


%% bv = bv.updateroi(analysis.stackch1,'polygon') update roi 
eps = roi(analysis2.stackch2,'polygon'); % eps: extraparenchymal space
eps = eps.modifyroi(analysis2.stackch1);
eps.stacks.ch1 = eps.addstack(analysis2.stackch1);
eps = eps.radonthresholding('ch1');
eps.stacks.ch2 = eps.addstack(analysis2.stackch2);
eps = eps.radonthresholding('ch2');


%%
%% bv = bv.updateroi(analysis.stackch1,'polygon') update roi 
eps2 = roi(analysis2.stackch2,'polygon'); % eps: extraparenchymal space
eps2 = eps2.modifyroi(analysis2.stackch1);
%%
eps2.stacks.ch1 = eps2.addstack(analysis2.stackch1);
eps2 = eps2.radonthresholding('ch1');
eps2.stacks.ch2 = eps2.addstack(analysis2.stackch2);
eps2 = eps2.radonthresholding('ch2');
%%
eps2.showstack('ch1')
eps2.showstack('tirs_ch1')
eps2.showstack('irtd_ch1')
%%
eps2.showroi

%%



tmp.stackch1 = analysis2.stackch1;
tmp.stackch2 = analysis2.stackch2;

%% Match the result of irdaon to the image stack
tmp.stackch1 = analysis2.stackch1;
% tmp.stackch1(tmp.minyb:tmp.maxyb, tmp.minxb:tmp.maxxb,:) = tmp.stackch1(tmp.minyb:tmp.maxyb, tmp.minxb:tmp.maxxb,:).*uint16(tmp.irtdthr);
tmp.stackch2 = analysis2.stackch2;
% tmp.stackch2(tmp.minyb:tmp.maxyb, tmp.minxb:tmp.maxxb,:) = tmp.stackch2(tmp.minyb:tmp.maxyb, tmp.minxb:tmp.maxxb,:).*uint16(tmp.irtdthr);

% Calculate half-size of the image stack (rounded up)
tmp.irtdthr = eps2.stacks.irtd_ch1<0.05;


%%
tmp.szirtd = ceil(size(eps2.stacks.irtd_ch1, 1) / 2);

% Calculate the midpoint of the ROI in x-direction
tmp.xhalf = (max(eps2.vertices(:, 1)) + min(eps2.vertices(:, 1)));
tmp.xhalf = ceil(tmp.xhalf / 2);
tmp.maxxb = tmp.xhalf + tmp.szirtd;
tmp.minxb = tmp.xhalf - tmp.szirtd + 2;

% Calculate the midpoint of the ROI in y-direction
tmp.yhalf = (max(eps2.vertices(:, 2)) + min(eps2.vertices(:, 2)));
tmp.yhalf = ceil(tmp.yhalf / 2);
tmp.maxyb = tmp.yhalf + tmp.szirtd ;
tmp.minyb = tmp.yhalf - tmp.szirtd + 2;

% Ensure boundaries are within image dimensions
[imgH, imgW] = size(tmp.stackch1(:, :, 1)); % Image dimensions
tmp.minxb = max(tmp.minxb, 1);
tmp.maxxb = min(tmp.maxxb, imgW);
tmp.minyb = max(tmp.minyb, 1);
tmp.maxyb = min(tmp.maxyb, imgH);

figure()
sliceViewer(tmp.stackch1)

%% convert original stack to rgb stack to apply mask

figure()
sliceViewer(tmp.stackch2)
tmp.rgb = zeros([size(tmp.stackch1),3],'uint8');
tmp.rgb(:,:,:,1) =  mat2gray(tmp.stackch1)*255;
tmp.rgb(:,:,:,2) = mat2gray(tmp.stackch2)*255;
tmp.rgb(tmp.minyb:tmp.maxyb, tmp.minxb:tmp.maxxb,:,3) = mat2gray(uint16(~tmp.irtdthr))*255;
%%
figure("Name",'rgb')
sliceViewer(tmp.rgb)
%%
io_savetiff(tmp.rgb,fullfile(analysis2.info.analysis_savepath,'rgb2.tiff'),tmp.save_resolution)
%%
%%
[tmp.xgrid,tmp.ygrid] = (ndgrid(1:1:size(eps2.stacks.tirs_ch1,1),1:1:size(eps2.stacks.tirs_ch1,2)));
tmp.row_coordinate = eps2.stacks.tirs_ch1.*tmp.xgrid;
tmp.radon2pi = squeeze(min(tmp.row_coordinate + (tmp.row_coordinate == 0) * 9999, [], 1));
tmp.radonpi = squeeze(max(tmp.row_coordinate,[],1));
tmp.radonmax = squeeze(sum(tmp.row_coordinate,1))-tmp.radonpi-tmp.radon2pi;

fwhm.theta = tmp.radonpi-tmp.radon2pi;

%%

tmp.fps = str2double(analysis2.info.savefps);
fwhm.taxis = linspace(0,size(fwhm.theta,2)/tmp.fps,size(fwhm.theta,2));

%%
figure('Name','diameter at theta')
imagesc(fwhm.taxis, 1:size(fwhm.theta, 1), fwhm.theta.*str2double(analysis2.info.objpix(1:end-2)));
xlabel('Time (s)');
ylabel('Theta (deg.)');
cb = colorbar; % Optional, to show the color scale

% Move the colorbar label to the top of the colorbar
title('Diameter at Theta');
%%
%%
fwhm.normtheta = fwhm.theta./mean(fwhm.theta,2);
fwhm.normstd = std(fwhm.theta,[],2);
%%
[tmp.maxstd,fwhm.mxlocstd] = max(fwhm.normstd);
[tmp.minstd,fwhm.mnlocstd] = min(fwhm.normstd);

%%
figure('Name','Normalized diameter at theta')
imagesc(fwhm.taxis, 1:size(fwhm.normtheta, 1), fwhm.normtheta-1);
cb = colorbar; % Optional, to show the color scale
hold on; % Keep the plot for overlaying the line
ylabel('angle (deg.)', 'Color', 'white', 'FontSize', 16);       % Label with font size
xlabel('time (sec)', 'Color', 'white', 'FontSize', 16);       % Label with font size

y_position = fwhm.mnlocstd; % Specify the y-coordinate of the horizontal line
yline(y_position, 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--'); % Add the horizontal line
y_position = fwhm.mxlocstd; % Specify the y-coordinate of the horizontal line
yline(y_position, 'LineWidth', 1, 'Color', 'magenta', 'LineStyle', '--'); % Add the horizontal line

x_position = 1060; % Specify the y-coordinate of the horizontal line
xline(x_position, 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--'); % Add the horizontal line
x_position = 1150; % Specify the y-coordinate of the horizontal line
xline(x_position, 'LineWidth', 1, 'Color', 'magenta', 'LineStyle', '--'); % Add the horizontal line
ax = gca;
ax.Color = 'black';           % Set background color of the axes
ax.XColor = 'white';          % Set x-axis color
ax.YColor = 'white';   
hold off;
%xlim([1060 1150])


%% standard deviation plot
% Create the figure
figure('Name', 'Inverted Plot with Custom Font Size');

% Plot the data
plot(fwhm.normstd, 'LineWidth', 1.5, 'Color', 'white'); % Make the plot line white

% Set axes properties
ax = gca;
ax.Color = 'black';           % Set background color of the axes
ax.XColor = 'white';          % Set x-axis color
ax.YColor = 'white';          % Set y-axis color
ax.LineWidth = 1.5;           % Optionally, increase the line thickness of axes
ax.FontSize = 14;             % Set the font size of the axes

% Add labels and title with custom font size
xlabel('angle (deg.)', 'Color', 'white', 'FontSize', 16);       % Label with font size
ylabel('Normalized Std', 'Color', 'white', 'FontSize', 16); % Label with font size
title('Normalized Std over Time', 'Color', 'white', 'FontSize', 18); % Title with font size

% Set the figure background
set(gcf, 'Color', 'black'); % Set the figure background color to black
xlim([0 180])
%%

% Circular indexing
tmp.maxindices = mod((fwhm.mxlocstd - 15:fwhm.mxlocstd + 15) - 1, size(fwhm.theta,1)) + 1;
fwhm.max30dstdmean = mean(fwhm.normtheta(tmp.maxindices,:),1);
%%
figure()
plot(medfilt1(fwhm.max30dstdmean-1,5))
%%
fwhm.filtered30stdmax = medfilt1(fwhm.max30dstdmean-1,50);


%% ball signal


% Define minimum stop duration (5 seconds)
tmp.min_stop_duration = 10;

% Extract time and velocity
time = analysis2.analog.ball_table{:, 1};

abs_velocity = abs(analysis2.analog.ball_table{:, 2});
abs_velocity = medfilt1(abs_velocity,1000);
abs_velocity(abs_velocity < 0.005) = 0;

figure('Name','absvelocity','NumberTitle','off')
plot(abs_velocity)

%% Plot locomotin trigger point
% Threshold to detect movement
tmp.logic_move = abs_velocity > 0.02;

% Identify continuous "stops" (gaps in movement)
[labelled_regions, num_regions] = bwlabel(~tmp.logic_move);  % Inverse to find stops

% Initialize filtered logic array
filtered_logic_move = tmp.logic_move;

% Process each stop region
for region_idx = 1:num_regions
    % Get indices for the current stop region
    stop_indices = find(labelled_regions == region_idx);
    
    % Calculate the duration of the stop
    tmp.stop_duration = time(stop_indices(end)) - time(stop_indices(1));
    
    % If stop duration is less than 5 seconds, replace it with running (set to 1)
    if tmp.stop_duration < tmp.min_stop_duration
        filtered_logic_move(stop_indices) = 1;
    end
end
trigger.locoffset = 5000;
trigger.transition = find(diff(filtered_logic_move)==1)-trigger.locoffset;

% Plot original and filtered signals
figure;
% subplot(4, 1, 1);
% plot(time, tmp.logic_move, 'r');
% title('Thresholded Movement');
% xlabel('Time (s)');
% ylabel('Movement (0 or 1)');
% xlim([0 1800]);
% 
% subplot(4, 1, 2);
% plot(time, filtered_logic_move, 'b');
% title('Continuous locomotion');
% xlabel('Time (s)');
% ylabel('Movement (0 or 1)');
% xlim([0 1800]);
% hold on
% plot(time(transition) + 1, 0, 'rx', 'MarkerSize', 8,'LineWidth', 2);  % Mark transitions

subplot(2,1,1)
plot(time,analysis2.analog.ball_table{:, 2});
hold on 
plot(time(transition) + 1, 0, 'rx', 'MarkerSize', 8,'LineWidth', 2);  % Mark transitions
plot(time(transition(trigger.thrloc)) + 1, 0, 'gx', 'MarkerSize', 8,'LineWidth', 2);  % Mark transitions
xlim([0 1800]);
title('Raw locomotion');

subplot(2,1,2)
plot(fwhm.taxis,fwhm.filtered30stdmax)
hold on
plot(fwhm.taxis(ceil(time(transition).*tmp.fps)), fwhm.filtered30stdmax(ceil(time(transition).*tmp.fps)), 'rx', 'MarkerSize', 10,'LineWidth', 2);  % Mark transitions
plot(fwhm.taxis(ceil(time(transition(trigger.thrloc)).*tmp.fps)), fwhm.filtered30stdmax(ceil(time(transition(trigger.thrloc)).*tmp.fps)), 'gx', 'MarkerSize', 10,'LineWidth', 2);  % Mark transitions

xlim([0 1800])

%% Locomotion triggred average
trigger.time = ceil(time(transition).*tmp.fps);

% eps.stacks.ch1(:,:,)
% analysis.stackch1
% 
% triggered_pvs = zeros([length(locpik_area)-2,501]);
trigger.lococh1 = zeros([size(analysis2.stackch1,1),size(analysis2.stackch1,2),ceil(10*tmp.fps),length(trigger.time)]);
trigger.locodiam = zeros(ceil(10*tmp.fps),length(trigger.time));
for idx = 1:length(trigger.time)
    disp(idx)
    trigger.lococh1(:,:,:,idx) = analysis2.stackch1(:,:,trigger.time(idx):trigger.time(idx)+ceil(10*tmp.fps)-1);
    trigger.lococh2(:,:,:,idx) = analysis2.stackch2(:,:,trigger.time(idx):trigger.time(idx)+ceil(10*tmp.fps)-1);
    trigger.locodiam(:,idx) = fwhm.filtered30stdmax(trigger.time(idx):trigger.time(idx)+ceil(10*tmp.fps)-1);
end

trigger.analogfreq = str2double(analysis2.info.analogfreq(1:end-2));
trigger.locotaxis = linspace(-trigger.locoffset/trigger.analogfreq,ceil(10*tmp.fps)/tmp.fps,1+ceil(10*tmp.fps)-trigger.locoffset/trigger.analogfreq);
%% Locomotion triggered average vasodilation
figure();
hold on;
% Plot the averaged line as a thick black line
plot(trigger.locotaxis,mean(trigger.locodiam, 2), 'k', 'LineWidth', 2); % Thick black line for the mean
tmp.fps
% Plot individual lines with low alpha (faded appearance)
for idx = 1:length(trigger.time)
    plot(trigger.locotaxis,trigger.locodiam(:, idx), 'Color', [0.8, 0.8, 0.8, 0.7]); % Faded gray lines with low alpha
end
xlim([0 10]);
ylim([-0.07 0.13]);
% Add labels and title for clarity (optional)
xlabel('Time (sec)');
ylabel('Norm. diameter');
title('Loco Diameter with Averaged Line');
grid off; % Add grid for better readability
%% filtering if locomotion triggered vessel diameter constrict after the trigger

tmp.lochalf = median(trigger.locodiam(1:ceil(end/2),:),1);
tmp.lochalf2 = median(trigger.locodiam(ceil(end/2):end,:),1);
tmp.locdiff = tmp.lochalf2-tmp.lochalf;
trigger.thrloc = find(tmp.locdiff > 0.05); 
trigger.locodiamfiltered = trigger.locodiam(:,trigger.thrloc);


%% plot filtered Locomotion triggered average vasodilation
figure();
hold on;
% Plot the averaged line as a thick black line
plot(trigger.locotaxis,mean(trigger.locodiamfiltered, 2), 'k', 'LineWidth', 2); % Thick black line for the mean
tmp.fps
% Plot individual lines with low alpha (faded appearance)
for idx = 1:length(trigger.thrloc)
    plot(trigger.locotaxis,trigger.locodiamfiltered(:, idx), 'Color', [0.8, 0.8, 0.8, 0.7]); % Faded gray lines with low alpha
end
xlim([0 10]);
ylim([-0.07 0.13]);
% Add labels and title for clarity (optional)
xlabel('Time (sec)');
ylabel('Norm. diameter');
title('> 5% Loco Diameter with Averaged Line');
grid off; % Add grid for better readability
% 
% figure('Name','Dilatation triggered average')
% plot(mean(triggered_pvs,1))
% xlim([0 500])

%%

trigger.locfilteredavgch1 = mean(trigger.lococh1(:,:,:,trigger.thrloc),4);
trigger.locfilteredavgch2 = mean(trigger.lococh2(:,:,:,trigger.thrloc),4);


%%
io_savetiff(trigger.locfilteredavgch1,fullfile(analysis2.info.analysis_savepath,'filtlocavgch1.tiff'),tmp.save_resolution)
%%
io_savetiff(trigger.locfilteredavgch2,fullfile(analysis2.info.analysis_savepath,'filtlocavgch2.tiff'),tmp.save_resolution)

%%
figure()
sliceViewer(ans)



%%

fwhm.mean = mean(fwhm.theta,2);
%%

[vertices,~] = roi_rectangle_polygon(eps2.stacks.irtd_ch1,'polygon');

%%
[~,maxloc] = max(fwhm.std);

%%
[~,minloc] = min(fwhm.std);
minindices = mod((minloc - 15:minloc + 15) - 1, size(fwhm.theta,1)) + 1;
fwhm.min10dstdmean = mean(fwhm.theta(minindices,:),1);
%%
figure()
plot(fwhm.min10dstdmean)

%%
fwhm.median = median(fwhm.theta,2);
%%
irtdsli = gpuArray(zeros(size(eps2.stacks.tirs_ch1,1),size(eps2.stacks.tirs_ch1,1),size(eps2.stacks.tirs_ch1,3)));
gputirs = gpuArray(eps2.stacks.tirs_ch1);
tic
for i = 1:size(gputirs,3) 
irtdsli(:,:,i) = iradon(gputirs(:,:,i),(1:1:180),'linear','None',1,size(eps2.stacks.tirs_ch1,1));
end
%%
toc
irtdsli = gather(irtdsli);
util_checkstack(irtdsli)
%%

figure()
%plot(prctile(fwhm.theta,50,1))
plot(fwhm.min10dstdmean)

%%
%%
figure()
imshow(row_coordinate(:,:,1))
%%
sinogram = false(size(eps2.stacks.tirs_ch1,1),size(eps2.stacks.tirs_ch1,2));
%%
%%
util_checkstack(eps2.stacks.tirs)
%%
sinogram = false(size(eps2.stacks.tirs_ch1,1),size(eps2.stacks.tirs_ch1,2));
radius1 = 20;
angle1 = 31;
sinogram(1:round(size(sinogram,1)/2)+radius1+5,angle1:angle1) = 1;
irtdsli = iradon(sinogram,(1:1:180),'pchip','None',1,size(eps2.stacks.tirs_ch1,1));
% fft2(irtdsli)
% fft(sinogram(:,angle1))

%irtdsli = iradon(eps.stacks.tirs_ch1(:,:,14600),(1:1:180),'linear','Hamming',1,size(eps.stacks.tirs_ch1,1));
figure()
subplot(1,2,1)
imshow(sinogram)
subplot(1,2,2)
imshow(irtdsli*100)
colormap(gca, 'Parula'); % Apply the Inferno colormap to the current axes

%%
% Step 1: Create the sinogram with a single projection line
sinogram = false(size(eps2.stacks.tirs_ch1,1), size(eps2.stacks.tirs_ch1,2));
radius1 = 20;
angle1 = 31; % Angle of the projection line
sinogram(round(size(sinogram,1)/2)+radius1:round(size(sinogram,1)/2)+radius1+5, angle1) = 1;

% Step 2: Inverse Radon transform (iradon)
irtdsli = iradon(sinogram, (1:1:180), 'linear', 'None', 1, size(eps2.stacks.tirs_ch1, 1));

% Step 3: Convert irtdsli to polar coordinates
[m, n] = size(irtdsli);
[X, Y] = meshgrid(1:n, 1:m);
center = [ceil(n/2), ceil(m/2)];
R = sqrt((X-center(1)).^2 + (Y-center(2)).^2);
Theta = atan2(Y-center(2), X-center(1));
r_max = max(R(:));
num_r = m; % Radial resolution
num_theta = n; % Angular resolution
r = linspace(0, r_max, num_r);
theta = linspace(-pi, pi, num_theta);
[R_polar, Theta_polar] = meshgrid(r, theta);
X_polar = R_polar .* cos(Theta_polar) + center(1);
Y_polar = R_polar .* sin(Theta_polar) + center(2);
polar_irtdsli = interp2(X, Y, irtdsli, X_polar, Y_polar, 'linear', 0);

% Step 4: Fourier transform comparisons
fft2_cartesian = fft2(irtdsli); % FFT2 of Cartesian reconstructed image
fft2_polar = fft2(polar_irtdsli); % FFT2 of polar reconstructed image
sinogram_fft = fft(sinogram(:, angle1)); % 1D FFT of the sinogram column

% Step 5: Display results
figure;

% Display Cartesian reconstructed image (irtdsli)
subplot(2, 3, 1);
imshow(irtdsli, []);
title('Cartesian irtdsli (Reconstructed Image)');

% Display Polar converted irtdsli
subplot(2, 3, 2);
imshow(polar_irtdsli, []);
title('Polar irtdsli');

% Display FFT2 of Cartesian irtdsli
subplot(2, 3, 3);
imshow(log(abs(fftshift(fft2_cartesian)) + 1), []);
title('FFT2 of Cartesian irtdsli');

% Display FFT2 of Polar irtdsli
subplot(2, 3, 4);
imshow(log(abs(fftshift(fft2_polar)) + 1), []);
title('FFT2 of Polar irtdsli');

% Display 1D FFT of sinogram projection
subplot(2, 3, 5);
plot(abs(sinogram_fft));
title('1D FFT of sinogram (angle 31)');

% Fourier Slice Theorem visualization
subplot(2, 3, 6);
imshow(log(abs(fftshift(sinogram_fft)) + 1), []);
title('Fourier Slice Theorem (1D FFT in 2D)');



%%

%%

figure()
imagesc(mean(eps2.stacks.tirs,3))


%%

%%
sinogram = false(size(eps2.stacks.tirs_ch1,1),size(eps2.stacks.tirs_ch1,2));
radius1 = 20;
radius2 = 10;
angle1 = 31;
angle2 = 30;
dif = 30;
sinogram(round(size(sinogram,1)/2)+radius1:round(size(sinogram,1)/2)+radius1,angle1-dif:angle1-dif) = 1;
sinogram(round(size(sinogram,1)/2)+radius1:round(size(sinogram,1)/2)+radius1,angle1:angle1) = 1;
sinogram(round(size(sinogram,1)/2)+radius1:round(size(sinogram,1)/2)+radius1,angle1+dif:angle1+dif) = 1;

%sinogram(round(size(sinogram,1)/2)+radius1:round(size(sinogram,1)/2)+radius1,angle2:angle2) = 1;

% sinogram(round(size(sinogram,1)/2)-radius1:round(size(sinogram,1)/2)-radius1,angle1-24:angle1-15) = 1;

% sinogram(round(size(sinogram,1)/2)-radius1:round(size(sinogram,1)/2)-radius1,angle2:angle2) = 1;

irtdsli = iradon(sinogram,(1:1:180),'linear','None',1,size(eps2.stacks.tirs_ch1,1));

%irtdsli = iradon(eps.stacks.tirs_ch1(:,:,14600),(1:1:180),'linear','Hamming',1,size(eps.stacks.tirs_ch1,1));
figure()
subplot(1,2,1)
imshow(sinogram)
subplot(1,2,2)
imshow(irtdsli*50)
colormap(gca, 'Parula'); % Apply the Inferno colormap to the current axes

%%

vis_img = masks{1} + masks2{5} + masks3{10};

figure()
subplot(1,2,1)
imshow(vis_img)
subplot(1,2,2)
imagesc(radon(vis_img,1:1:360))
%%

%%
iradonsli = iradon(sinogram,(1:1:180),'linear','Hamming',1,size(eps2.stacks.tirs,1));
figure()
imagesc(irdadon)

%%
analog = analyze_readanalog(info.analyzefolder);
zstack = struct();
zstack.ch1 = analyze_readtiff(info.analyzefolder,'*ch1.tif');
zstack.ch2 = analyze_readtiff(info.analyzefolder,'*ch2.tif');
%% make x axis
fps = str2double(info.savefps);
fname.zstack = fieldnames(zstack);
taxis = linspace(0,size(zstack.(fname.zstack{1}),3)/fps,size(zstack.(fname.zstack{1}),3));



%% Median and gaussian filter 2D stack
filteredStack = gpuArray(zstack.ch1);
for i = 1:size(filteredStack, 3)
    sli = medfilt2(filteredStack(:, :, i), [3 3]);  % Apply 3x3 median filter (adjust size as needed)
    filteredStack(:, :, i) = imgaussfilt(sli,1);
end
filteredStack = gather(filteredStack);
filteredStack = pre_thresholding(filteredStack);

figure('Name','filtered STack')
sliceViewer(filteredStack)

%%
% vessel area selection change to polygon selection might beneficial once
% intracellular signal start to collected
[roi.bv_vertices, ~] = roi_rectangle_polygon(filteredStack,'polygon'); %2.1
roi.bv_vertices = round(roi.bv_vertices); 
roi.bv_stack = roi_applyvertices(filteredStack,roi.bv_vertices);
%% vessel lumen area calculation
tirs = analyze_radon(bv_stack);
%%
figure()
sliceViewer(tirs.irtd_norm)
%%
edge_image = edge(tirs.irtd_norm(:,:,1), 'Canny');          % Detect edges using the Canny method
figure()
imagesc(edge_image)
%%
[H, theta, rho] = hough(edge_image);
figure()
imshow(imadjust(rescale(H)), 'XData', theta, 'YData', rho, ...
       'InitialMagnification', 'fit');
xlabel('\theta (degrees)');
ylabel('\rho (pixels)');
title('Hough Transform (Parameter Space)');
%%
% Detect circles using Hough Circle Transform
[centers, radii] = imfindcircles(edge_image, [10 50], 'ObjectPolarity', 'bright', 'Sensitivity', 0.9);

% Visualize the detected circles
figure;
imshow(edge_image);
hold on;
viscircles(centers, radii, 'Color', 'b'); % Overlay detected circles
title('Detected Circles');

% Create a mask for the detected circle
circle_mask = false(size(edge_image));
[x_grid, y_grid] = meshgrid(1:size(edge_image, 2), 1:size(edge_image, 1));
for i = 1:size(centers, 1) % In case there are multiple circles
    distance_from_center = sqrt((x_grid - centers(i, 1)).^2 + (y_grid - centers(i, 2)).^2);
    circle_mask = circle_mask | (distance_from_center <= radii(i));
end

% Extract only the circle
circle_only_image = edge_image & circle_mask;

% Display results
figure;

subplot(1, 2, 1);
imshow(edge_image);
title('Original Edge Image');

subplot(1, 2, 2);
imshow(circle_only_image);
title('Circle Only');
%%
% Detect peaks in the Hough Transform
peaks = houghpeaks(H, 5, 'Threshold', ceil(0.3 * max(H(:)))); % Adjust parameters as needed

% Extract lines from Hough Transform
lines = houghlines(edge_image, theta, rho, peaks, 'FillGap', 20, 'MinLength', 30);

% Visualize the results
figure;


% Plot detected lines
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];     % Line endpoints
    plot(xy(:, 1), xy(:, 2), 'LineWidth', 2, 'Color', 'red'); % Plot line
    % Plot line endpoints
    plot(xy(1, 1), xy(1, 2), 'x', 'LineWidth', 2, 'Color', 'yellow');
    plot(xy(2, 1), xy(2, 2), 'x', 'LineWidth', 2, 'Color', 'green');
end

title('Detected Lines Using Hough Transform');
hold off;





%% Add 
fwhm_diam = squeeze(sum(tirs.tirs,1));

%%

norm_fwhmdiam = fwhm_diam-min(fwhm_diam,[],2);
norm_fwhmdiam = norm_fwhmdiam./max(fwhm_diam,[],2);
%% unevenness of dialation  
std_diam = medfilt1(squeeze(std(norm_fwhmdiam,[],1)),100);
mean_diam = medfilt1(squeeze(mean(norm_fwhmdiam,1)),100);
%%
norm_stddiam = std_diam-min(std_diam);
norm_stddiam = norm_stddiam/max(norm_stddiam);
norm_meandiam = mean_diam-min(mean_diam);
norm_meandiam = norm_meandiam/max(norm_meandiam);

%%
figure()
plot(norm_meandiam)
hold on;

%%
fill(taxis)

%%
figure()
plot(norm_stddiam)
hold on


%%
figure()
plot(medfilt1(squeeze(mean(norm_fwhmdiam,1)),100))

%%
figure()
imagesc(norm_fwhmdiam)
%%
median(fwhm_diam,2)

%% area calculation
disp('Area calculation and contour start')
tirs.area = zeros([size(tirs.irtd_norm,3),1]);
tirs.radoncontours=cell(size(tirs.irtd_norm,3),1); %countour lines obtained of the vessel lumen using TiRS method
%%
f=50;
irtd_threshold = 0.2;
[cc,l]=bwboundaries(tirs.irtd_norm(:,:,f)>irtd_threshold*max(tirs.irtd_norm(:,:,f),[],'all'));
numPixels = cellfun(@length,cc);
[~,idx] = max(numPixels);
area_filled=regionprops(l,'FilledArea','Image','FilledImage');
tirs.area(f)=length(find(area_filled(idx).FilledImage));
tirs.radoncontours{f} = contour(tirs.irtd_norm(:,:,f),[irtd_threshold irtd_threshold]*max(tirs.irtd_norm(:,:,f),[],'all'),'r', 'LineWidth', 2);
subplot(2,1,1)
imshow(l);
subplot(2,1,2)
%%

% Fill the Z-stack
for i = 1:num_frames
    % Get vertices for the current frame
    vertices = vertices_cell{i};
    
    % Create a binary mask for the current frame
    mask = poly2mask(vertices(1, :), vertices(2, :), image_height, image_width);
    
    % Store the mask in the Z-stack
    zstack(:, :, i) = mask;
end

% Display a sample frame
figure;
imshow(zstack(:, :, 1)); % Show the first frame
title('First Frame of Z-stack');

%% a
figure()

for f = 1:size(tirs.irtd_norm,3)
    [cc,l]=bwboundaries(tirs.irtd_norm(:,:,f)>irtd_threshold*max(tirs.irtd_norm(:,:,f),[],'all'));
    numPixels = cellfun(@length,cc);
    [~,idx] = max(numPixels);
    area_filled=regionprops(l,'FilledArea','Image','FilledImage');
    tirs.area(f)=length(find(area_filled(idx).FilledImage));
    tirs.radoncontours{f} = contour(tirs.irtd_norm(:,:,f),[irtd_threshold irtd_threshold]*max(tirs.irtd_norm(:,:,f),[],'all'),'r', 'LineWidth', 2);
    imshow(l);
end
disp('Area calculation and contour end')




%%
figure()
plot(analog.ball_table{:,1},abs_velocity);


%%




%%

% Initialize filtered logic array
filtered_logic_move = logic_move;

% Process each region
for region_idx = 1:num_regions
    % Get indices for the current region
    region_indices = find(labelled_regions == region_idx);
    
    % Calculate the duration of the region
    region_duration = time(region_indices(end)) - time(region_indices(1));
    
    % If duration is less than 5 seconds, remove it
    if region_duration < min_duration
        filtered_logic_move(region_indices) = 1;
    end
end



%%

abs_velocity = abs(analog.ball_table{:,2});
%%
logic_move = abs_velocity;
logic_move(logic_move<0.01)=0;
logic_move(logic_move>0.01)=1;
%%

plot(analog.ball_table{:,1},logic_move)





%%

vesselspectrum = analog_ecogspectrum(fps,filtered_area);

%%
figure()
pspectrum(filtered_area-mean(filtered_area),fps)
%%
%%
figure()
sliceViewer(tirs.irtd_norm)


%% Peak estimation 
filtered_area = medfilt1(tirs.area,10);
filtered_area = lowpass(filtered_area,2,12.4);

peak_area = filtered_area;
base_line = median(filtered_area);
peak_area(peak_area<base_line*1.3)=0;
%%
figure('Name','area plot')
plot(peak_area.^0.5)
%%
[pik_area,locpik_area] = findpeaks(peak_area,1,'MinPeakDistance',10,'MinPeakProminence',base_line);

figure('Name','Vessel area')
plot(filtered_area.^0.5)
% hold on
% plot(locpik_area,pik_area,'o')


%%
[pvs_vertices, ~] = roi_rectangle_polygon(pre_groupaverage(zstack.ch2,50),'polygon');

pvs_stack = roi_applyvertices(zstack.ch2,pvs_vertices);

pvs_signal = squeeze(mean(pvs_stack,[1,2],'omitnan'));


%%
figure('Name','Individual PVS trace')
for idx = 1:20
plot(pvs_signal(locpik_area(idx)-250:locpik_area(idx)+250))
hold on
end
xlim([0 500])

%%
triggered_pvs = zeros([length(locpik_area)-2,501]);
for idx = 3:10 %1:length(locpik_area)-2
triggered_pvs(idx,:) = pvs_signal(locpik_area(idx+2)-250:locpik_area(idx+2)+250);
end

figure('Name','Dilatation triggered average')
plot(mean(triggered_pvs,1))
xlim([0 500])




%% ROI signal1
[pvs_vertices, ~] = roi_rectangle_polygon(pre_groupaverage(zstack.ch2,50),'polygon');

pvs_stack = roi_applyvertices(zstack.ch2,pvs_vertices);

pvs_signal = squeeze(mean(pvs_stack,[1,2],'omitnan'));
figure('Name','Individual PVS trace')
for idx = 1:20
plot(pvs_signal(locpik_area(idx)-250:locpik_area(idx)+250))
hold on
end
xlim([0 500])

%
triggered_pvs = zeros([length(locpik_area)-2,501]);
for idx = 3:10 %1:length(locpik_area)-2
triggered_pvs(idx,:) = pvs_signal(locpik_area(idx+2)-250:locpik_area(idx+2)+250);
end

figure('Name','Dilatation triggered average')
plot(mean(triggered_pvs,1))
xlim([0 500])
%% ROI signal2
[pvs_vertices, ~] = roi_rectangle_polygon(pre_groupaverage(zstack.ch2,50),'polygon');

pvs_stack = roi_applyvertices(zstack.ch2,pvs_vertices);

pvs_signal = squeeze(mean(pvs_stack,[1,2],'omitnan'));

figure()
plot(pvs_signal)


figure('Name','Individual ISF trace')
for idx = 3:length(locpik_area)
plot(pvs_signal(locpik_area(idx)-250:locpik_area(idx)+250))
hold on
end
xlim([0 500])



triggered_pvs = zeros([length(locpik_area)-2,501]);
for idx = 1:length(locpik_area)-2
triggered_pvs(idx,:) = pvs_signal(locpik_area(idx+2)-250:locpik_area(idx+2)+250);
end

figure('Name','Dilatation triggered ISF average')
plot(mean(triggered_pvs,1))
xlim([0 500])
%%
triggered_zstack= zeros(size(zstack.ch2,1),size(zstack.ch2,2),101,length(locpik_area)-10);

for idx = 1:length(locpik_area)-10
    
    triggered_zstack(:,:,:,idx) = zstack.ch2(:,:,locpik_area(idx+2)-50:locpik_area(idx+2)+50);
end
%%
roi_rectangle_polygon(mean(triggered_zstack,4))

%%
triggered_zstackbackground= zeros(size(zstack.ch2,1),size(zstack.ch2,2),201,length(locpik_area));
exclusion = 3;
for idx = 1:length(locpik_area)-exclusion    
    triggered_zstackbackground(:,:,:,idx) = zstack.ch2(:,:,locpik_area(idx+exclusion)-100:locpik_area(idx+exclusion)+100);
end
triggered_avgimg = mean(triggered_zstackbackground,4);
figure()
sliceViewer(triggered_avgimg)
%% base line image

upper5 = find(filtered_area<base_line);
bottom5 = find(filtered_area>base_line*0.95);
baseframes = intersect(upper5,bottom5);
base_ch2 = zstack.ch2(:,:,baseframes);
base_ch1 = zstack.ch1(:,:,baseframes);
figure()
sliceViewer(pre_groupaverage(base_ch2,10))
%%
base_ch2img = mean(base_ch2,3)./max(base_ch2,[],'all')*5;
base_ch1img = mean(base_ch1,3)./max(base_ch1,[],'all')*5;

rgb_img = cat(3,base_ch1img,base_ch2img);
rgb_img = cat(3,rgb_img,zeros(size(base_ch1img)));

%%
figure()
imshow(rgb_img,[0 0.5])




%% Normalize
dff = (triggered_avgimg-base_ch2img)./base_ch2img;

%%
figure()
cmap = jet(256);

sliceViewer((triggered_avgimg-base_ch2img)./base_ch2img,"Colormap",cmap);
%%

mfac = 1;
cmap = jet(256);

figure('Name', 'during');
imshow(mean(dff(:, :, 99:104), 3) * mfac, [0 1]);  % Display with range [0 1]
colormap(cmap);  % Apply the colormap
colorbar;  % Optional: Add a colorbar for reference

%%


figure('Name','before')
imshow(mean(dff(:, :, 1:6), 3) * mfac, [0 1]);  % Display with range [0 1]
colormap(cmap);  % Apply the colormap
colorbar;  % Optional: Add a colorbar for reference


figure('Name','after')
imshow(mean(dff(:, :, 195:201), 3) * mfac, [0 1]);  % Display with range [0 1]
colormap(cmap);  % Apply the colormap
colorbar;  % Optional: Add a colorbar for reference








%% ROI signal2
[pvs_vertices, ~] = roi_rectangle_polygon(triggered_avgimg,'polygon');

roi_stack = roi_applyvertices(dff,pvs_vertices);

dff_local = squeeze(mean(roi_stack,[1,2],'omitnan'));

figure()
plot(dff_local)
xlim([0 201])

%%
figure('Name','Individual ISF trace')
for idx = 3:length(locpik_area)
plot(pvs_signal(locpik_area(idx)-250:locpik_area(idx)+250))
hold on
end
xlim([0 500])
%%


triggered_pvs = zeros([length(locpik_area)-2,501]);
for idx = 1:length(locpik_area)-2
triggered_pvs(idx,:) = pvs_signal(locpik_area(idx+2)-250:locpik_area(idx+2)+250);
end

figure('Name','Dilatation triggered ISF average')
plot(mean(triggered_pvs,1))
xlim([0 500])
%%

