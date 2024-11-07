clear, clc
%make mcsx obj, get general (info), two photon scanning microscope imaging (info_ tpsm), and imaging mode specific (info_mode)  
[info,  mobj] = io_initmdf();
%%

info_mode.analogcount   = str2double(mobj.ReadParameter('Analog Sample Count'));
info_mode.analogresolution   = mobj.ReadParameter('Analog Resolution');
info_mode.analogfreq    = mobj.ReadParameter('Analog Acquisition Frequency (Hz)');
info_mode.analogfreq    = str2double(info_mode.analogfreq(1:end -3));
info_mode.analogch4range = mobj.ReadParameter('Analog Ch 4 Input Range');
analog = double(mobj.ReadAnalog(1,info_mode.analogcount,0));

%%
downsample(analog,100)
%%
plot(ans)
%%
param.start         = 0;       % [sec]
param.duration      = 300;    % [sec]
param.pshift        = 1;        % [pixel]
param.medfilt       = [3 3 1]; % 3d median filter voxel size [Y X Z]
param.framestart  = round(1+param.start*info.fps); % fps*sec
param.totalframe = round(param.duration*info.fps); % fps*sec
param.frameend = param.totalframe + param.framestart;

if param.totalframe > info.fcount % if total frame to load exceed end of frame, load from start to the end
    param.frameend = info.fcount;
end    

zstack = io_readframes(mobj,1,[param.framestart, param.framestart+param.totalframe]); % read frame


% Preprcocessing (Padding removal -> post pixel shift correction -> Trim)
 
[tmp.xpadStart,tmp.xpadEnd] = pre_findpadding(zstack); % Find padded region caused by sinusoidal correction

% motion correction by dst registration
dst_stack = pre_groupaverage(zstack(:,tmp.xpadStart:tmp.xpadEnd,:), 10); % denoise by group averaging 
dst_stack = medfilt3(dst_stack,[3,3,5]); % denoise by 3d median filter
dst_stack = dst_stack - min(dst_stack,[],'all'); % min subtraction to make non negative matrix
drift_table = pre_estimatemotion(dst_stack); % using dft_registration, get drift table
%%
ampshiftx = max(pixelShift(4,:)) - min(pixelShift(4,:));
ampshifty = max(pixelShift(3,:)) - min(pixelShift(3,:));
motioncorrected_zstac = zeros([size(,)]);
for sli = drange(1:size(hold_stack,3))
    disp(sli)
end
max(pixelshift(3,:)
%%
[MScanData.notes.pixelShift(:,theFrame),~] = DftRegistration_2P(fftFirstFrame,fftRawFrame,1);





%%
figure(5)
sliceViewer(medstack,"DisplayRangeInteraction",'on')
%% Radon transform diameter analysis

maxframe = 10000000;
nframes=param.totalframe;% number of frames in file
the_angles=1:1:180;
rtd_threhold=.5;%. threshold in radon space, typically between 0.35 to 0.5
irtd_threhold=.2;%.threshold for converting back using inverse Radon transform, typically 0.2
area=-ones(min(nframes,maxframe),1); % area of vessel cross-section
radoncontours=cell(min(nframes,maxframe),1); %countour lines obtained of the vessel lumen using TiRS method

%displays the first frame of the .tif stack and asks the user to draw a box
%around the penetrating vessel

%% first 20frame averaged image 
%%
figure(77)
% display the first frame of the .tif stack so the user can draw a box around the vessel
imagesc(medstack(:,:,1))
colormap gray
axis image
title('draw a box around region to be processed')
thebox=drawrectangle();
input('type enter when done','s');
xy=round(thebox.Vertices);
raw_hold_image=double(medstack);
hold_stack=raw_hold_image(xy(1,2):xy(3,2), xy(1,1):xy(3,1),:);
%%
[path.tifName, path.tifPath] = uigetfile({'*.tif ; *.mdf'}); % select file by UI
thefile = [path.tifPath, path.tifName];
thefile_info=imfinfo(thefile);% the name of the file to be loaded
nframes=length(thefile_info);% number of frames in file
the_angles=1:1:180;
rtd_threhold=.5;%. threshold in radon space, typically between 0.35 to 0.5
irtd_threhold=.2;%.threshold for converting back using inverse Radon transform, typically 0.2
area=-ones(min(nframes,maxframe),1); % area of vessel cross-section
radoncontours=cell(min(nframes,maxframe),1); %countour lines obtained of the vessel lumen using TiRS method

%displays the first frame of the .tif stack and asks the user to draw a box
%around the penetrating vessel

figure(77)
% display the first frame of the .tif stack so the user can draw a box around the vessel
imagesc(imread(thefile,'tif', 'Index',1))
colormap gray
axis image
title('draw a box around region to be processed')
thebox=imrect(gca);
theinput=input('type enter when done','s');

api = iptgetapi(thebox);
mv_mpP(1).Vessel.box_position.xy=api.getPosition();
%%
f=1;
%read in the image in the tiff stack and converti it into a double for procesing
raw_hold_image=double(sum(imread(thefile,'tif', 'Index',f),3));
hold_image=raw_hold_image(round(mv_mpP(1).Vessel.box_position.xy(2):mv_mpP(1).Vessel.box_position.xy(2)+mv_mpP(1).Vessel.box_position.xy(4)),...
    round(mv_mpP(1).Vessel.box_position.xy(1):mv_mpP(1).Vessel.box_position.xy(1)+mv_mpP(1).Vessel.box_position.xy(3)));
radon_hold_image=radon(hold_image,the_angles);
%%
figure(4)
imshow(radon_hold_image,[0 100000])

%%
for f=1:min(nframes,maxframe)
    %read in the image in the tiff stack and converti it into a double for procesing
    raw_hold_image=double(sum(imread(thefile,'tif', 'Index',f),3));
    hold_image=raw_hold_image(round(mv_mpP(1).Vessel.box_position.xy(2):mv_mpP(1).Vessel.box_position.xy(2)+mv_mpP(1).Vessel.box_position.xy(4)),...
        round(mv_mpP(1).Vessel.box_position.xy(1):mv_mpP(1).Vessel.box_position.xy(1)+mv_mpP(1).Vessel.box_position.xy(3)));
    radon_hold_image=radon(hold_image,the_angles);
    for k=1:length(the_angles)
        % for each angle, rescale the image intensity to 0 to 1.
        radon_hold_image(:,k)=radon_hold_image(:,k)-min(radon_hold_image(:,k));
        radon_hold_image(:,k)=radon_hold_image(:,k)/max(radon_hold_image(:,k));
        [maxpoint(k),maxpointlocation(k)]=max(radon_hold_image(:,k));

        % find the threshold crossings in Radon Space
        [~,min_edge(k)]=max(find(radon_hold_image(1:maxpointlocation(k),k)<rtd_threhold));
        [~,max_edge(k)]=max(find(radon_hold_image(maxpointlocation(k)+1:end,k)>rtd_threhold));
        
        % set all other pixels to 0
        radon_hold_image(1:min_edge(k),k)=0;
        radon_hold_image((maxpointlocation(k)+max_edge(k)):end,k)=0;
    end
    %transform the image that has been threholded in Radon space back to
    %image space
    irtd_norm=iradon(double(radon_hold_image>rtd_threhold*max(radon_hold_image(:))),(the_angles),'linear','Hamming',size(radon_hold_image,2));
    
    %threshold the image
    [cc,l]=bwboundaries(irtd_norm>irtd_threhold*max(irtd_norm(:)));
    numPixels = cellfun(@length,cc);
    %find the largest contiguous group of suprathreshold pixels
    [~,idx] = max(numPixels);
    figure(44)
    
    subplot(221)
    hold off
    
    imagesc(hold_image,'XData',[1:size(hold_image,2)]+90-size(hold_image,2)/2,'YData',[1:size(hold_image,1)]+90-size(hold_image,1)/2)
    axis([min([1:size(hold_image,2)]+90-size(hold_image,2)/2) max([1:size(hold_image,2)]+90-size(hold_image,2)/2)...
        min([1:size(hold_image,1)]+90-size(hold_image,1)/2) max([1:size(hold_image,1)]+90-size(hold_image,1)/2)])
    axis equal
    
    axis xy
    hold on
    colormap gray
    
    title(['frame=' num2str(f)])
    subplot(222)
    imagesc(irtd_norm)
    axis image
    axis xy
    title('inverse-Radon transformed image')
    
    area_filled=regionprops(l,'FilledArea','Image','FilledImage');
    area(f)=length(find(area_filled(idx).FilledImage));
    subplot(221)
    radoncontours{f}=contour(irtd_norm(1:end,1:end),[irtd_threhold irtd_threhold]*max(irtd_norm(:)),'r', 'LineWidth', 2);
    pause(0.05)
end

%plot the area for each frame
subplot(2,1,2)
plot(1:f,area(1:f),'ro')
xlabel('frame number')
ylabel('area, pixels')
ylim([0 1.2*max(area(1:f))])