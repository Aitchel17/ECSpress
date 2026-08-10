function sliderseg(img)
% SLIDERSEG is a helper function to allow the user to preview the effect of
% various input parameters on the segmentation process in seg.m. See
% SegmentStack.m for a description of the segmentation process. Called by
% SegmentStack.m, and calls seg.m.
% img is a M X N X S stack of gray scale images
% Written by Kimberly Boster, kboster@ur.rochester.edu

    if isa(img,'uint16') || isa(img,'uint8')
        maxI=double(max(img(:)));
        img=double(img)/maxI;
    elseif isa(img,'logical')
        img=double(img);
        maxI=1;
    else
        maxI=1;
    end

    if size(img,3)>100 % if the slider bar gets too small, it throws an error. Therefore, if there are more than 100 images, downsample img so the slider doesn't get too skinny. 
        downsample_inc=ceil(size(img,3)/100);
        img=img(:,:,1:downsample_inc:size(img,3));
        disp(['img was downsampled for the slider, only keeping every ' num2str(downsample_inc) ' images'])
    end
    
    % control row locations
    rowfancythresh=.1;
    rowsmooth=0;
    rowsize=.05;
    rowthresh=.15;
    rowslice=.25;
    rowsat=.2;

    %initialize
    s=1;
    T=0.01;
    sat=1;
    satdisp=sat*maxI;
    Tdisp=T*sat*maxI;
    sizefilt=0;
    sizesmooth=0;
    fillcleanopt=0;
    autothreshopt=0;
    adaptthreshopt=0;
    dilateadapt=0;
    se=strel('disk',dilateadapt);
    disptype='overlay';
    
    figure,
    if length(size(img))>2
        sliceslider=superSlider(gcf, 'position',[0.15,rowslice,0.75,0.03],'numslides',1,'max',1,'min',1/size(img,3),'stepsize',1/size(img,3),'Callback',@sliderCallbackSlice);   
        sliceTextbox = uicontrol('Style','edit','units','normalized', 'Position',[0.92, rowslice,.08 .05],'String',num2str(s), 'Callback',@inputFromSliceTextbox);
    end
    threshslider=superSlider(gcf, 'position',[0.15,rowthresh,0.75,0.03],'numslides',1,'max',1,'min',0,'stepsize',.01,'Callback',@sliderCallbackThresh);  
    satslider=superSlider(gcf, 'position',[0.15,rowsat,0.75,0.03],'numslides',1,'max',1,'min',0,'stepsize',.01,'Callback',@sliderCallbackSat);
    sizefiltslider=superSlider(gcf, 'position',[0.15,rowsize,0.75,0.03],'numslides',1,'max',1,'min',0,'stepsize',.02,'Callback',@sliderCallbackSizeFilt); 
    sizesmoothslider=superSlider(gcf, 'position',[0.15,rowsmooth,0.5,0.03],'numslides',1,'max',1,'min',0,'stepsize',.05,'Callback',@sliderCallbackSizeSmooth); 
    dilateadaptslider=superSlider(gcf, 'position',[0.45,rowfancythresh,0.25,0.03],'numslides',1,'max',1,'min',0,'stepsize',.05,'Callback',@sliderCallbackDilateAdapt); 
    FillCleanOptBox = uicontrol('style','checkbox','units','normalized','position',[0.75,rowsmooth,0.25,0.03], 'string','fill holes, clean?','callback',@checkboxFillCleanOpt);
    autothreshOptBox = uicontrol('style','checkbox','units','normalized','position',[0.82,rowfancythresh,.2,0.03], 'string','auto thresh?','callback',@checkboxautothreshOpt);
    adaptthreshOptBox = uicontrol('style','checkbox','units','normalized','position',[0,rowfancythresh,.32,0.03], 'string','spatially varying thresh?','callback',@checkboxadaptthreshOpt);
    satTextbox = uicontrol('Style','edit','units','normalized', 'Position',[0.92, rowsat,.08 .05],'String',num2str(satdisp), 'Callback',@inputFromSatTextbox);
    threshTextbox = uicontrol('Style','edit','units','normalized', 'Position',[0.92, rowthresh,.08 .05],'String',Tdisp, 'Callback',@inputFromThreshTextbox);
    sizefiltTextbox = uicontrol('Style','edit','units','normalized', 'Position',[0.92, rowsize,.08 .05],'String',sizefilt, 'Callback',@inputFromSizefiltTextbox);
    dispmenu = uicontrol('Style','popupmenu','units','normalized', 'Position',[0,.95,.2 .05],'String',{'overlay','outline','image alone'}, 'Callback',@inputFromDispMenu);
    
   
    % update satslider location
    allSlides = get(satslider, 'Children');
    location = get(allSlides(1), 'Position');
    location(1)=sat*(1-2*location(3))+location(3); % scale since the scale bar has a finite width location(3) 
    set(allSlides(1),'Position',location) 
    
   
    updateimage()    
    disp('Choose slice, threshold, and filter size using slider. (Draging the bar does not work.)') %I'm not sure why dragging the bar doesn't work, since you can drag the bar on sliderbw_2thresh.

    function sliderCallbackSlice(~, ~)
        infoMatrix = get(sliceslider, 'UserData');
        s=round(infoMatrix(1,1)*size(img,3));
        updateimage()
        set(sliceTextbox,'String', num2str(s,'%0.5g'))
    end

    function sliderCallbackThresh(~, ~)
        infoMatrix = get(threshslider, 'UserData');
        T=infoMatrix(1,1);
        
        adaptthreshopt=0;
        set(adaptthreshOptBox,'value',0)
        autothreshopt=0;
        set(autothreshOptBox,'value',0)
        
        updateimage()      
    end

    function sliderCallbackSat(~, ~)
        infoMatrix = get(satslider, 'UserData');
        sat=infoMatrix(1,1);  
        satdisp=sat*maxI;
        set(satTextbox,'String', num2str(satdisp,'%0.5g'))
        
        updateimage()      
    end

    function inputFromSatTextbox(~,~)
        satdisp=str2double(get(satTextbox,'string')); 
        sat=satdisp/maxI;

        updateimage()  
        
        %update slider location
        allSlides = get(satslider, 'Children');
        location = get(allSlides(1), 'Position');
        location(1)=sat*(1-2*location(3))+location(3); % scale since the scale bar has a finite width location(3) 
        set(allSlides(1),'Position',location) 
    
    end 

    function inputFromSliceTextbox(~,~)
        s=round(str2double(get(sliceTextbox,'string'))); 
        
        updateimage()  
        
        %update slider location
        allSlides = get(sliceslider, 'Children');
        location = get(allSlides(1), 'Position');
        location(1)=(s/size(img,3))*(1-2*location(3))+location(3); % scale since the scale bar has a finite width location(3) 
        set(allSlides(1),'Position',location) 
    
    end 

    function inputFromSizefiltTextbox(~,~)
        sizefilt=round(str2double(get(sizefiltTextbox,'string'))); 
        
        updateimage()  
        
        %update slider location
        allSlides = get(sizefiltslider, 'Children');
        location = get(allSlides(1), 'Position');
        location(1)=sizefilt/50*(1-2*location(3))+location(3); % scale since the scale bar has a finite width location(3) 
        set(allSlides(1),'Position',location) 
    
    end 

    function inputFromThreshTextbox(~,~)
        Tdisp=get(threshTextbox,'string'); 
        T=str2double(Tdisp)/(maxI*sat);
        
        updateimage()  
        
        %update slider location
        allSlides = get(threshslider, 'Children');
        location = get(allSlides(1), 'Position');
        location(1)=T*(1-2*location(3))+location(3); % scale since the scale bar has a finite width location(3) 
        set(allSlides(1),'Position',location) 
    
    end 

    function sliderCallbackSizeFilt(~, ~)
        infoMatrix = get(sizefiltslider, 'UserData');
        sizefilt=50*infoMatrix(1,1);
        updateimage()
        set(sizefiltTextbox,'String', num2str(sizefilt,'%0.5g'))
    end

    function sliderCallbackDilateAdapt(~, ~)
        infoMatrix = get(dilateadaptslider, 'UserData');
        dilateadapt=round(20*infoMatrix(1,1));
        se=strel('disk',dilateadapt);
        updateimage()
    end

    function sliderCallbackSizeSmooth(~, ~)
        infoMatrix = get(sizesmoothslider, 'UserData');
        sizesmooth=round(20*infoMatrix(1,1));
        updateimage()
    end

    function checkboxFillCleanOpt(~, ~)
        fillcleanopt=get(FillCleanOptBox,'value');
        updateimage()
    end

    function checkboxautothreshOpt(~, ~)
        autothreshopt=get(autothreshOptBox,'value');
        if autothreshopt
            adaptthreshopt=0;
            set(adaptthreshOptBox,'value',0)
        else
            T=.01;
        end
        
        updateimage()
    end

    function checkboxadaptthreshOpt(~, ~)
        adaptthreshopt=get(adaptthreshOptBox,'value');
        if adaptthreshopt
            autothreshopt=0;
            set(autothreshOptBox,'value',0)
        else
            T=0.01;
        end
        updateimage()
    end
    function inputFromDispMenu(~, ~)
        val=get(dispmenu,'value');
        str=get(dispmenu,'string');
        disptype=str{val};
        updateimage()
    end

    function updateimage()
        if autothreshopt
            T=graythresh(img(:,:,s));
        end
        if adaptthreshopt
            T=adaptthresh(img(:,:,s)) + ~imdilate(imbinarize(img(:,:,s)),se);
            T(T>1)=1;
            Tdisp='adapt';
        else
            Tdisp=num2str(T*sat*maxI,'%0.5g');
        end
        set(threshTextbox,'String', Tdisp)
        
        imdisp=img(:,:,s)/sat;
        [imbw,BW_smooth]=seg(imdisp,T,sizefilt,sizesmooth,fillcleanopt);

        subplot('position',[0 .3 .6 .6]), 

        if strcmp(disptype,'overlay')
            imagesc(labeloverlay(imdisp,BW_smooth)), axis off % show original overlayed with segmentation 
            title('Final Overlaid on image')
        elseif strcmp(disptype,'outline')
            imagesc(imdisp, [0 1]), axis off, colormap('gray')
            BW_edge=bwperim(BW_smooth);
            [y,x]=ind2sub(size(BW_smooth),find(BW_edge));
            hold on, plot(x,y,'g.') % original + outline
            title('Outline of Segmentation')
        else
            imagesc(imdisp, [0 1]), colormap('gray'), axis off
            title('Original Image')
        end
     

        subplot('position',[.65 .65 .3 .3]),imshow(imbw),title('After threshold, filling, cleaning')
        subplot('position',[.65 .3 .3 .3]),imshow(BW_smooth),title('After smoothing & filtering')
        
        axtext=axes('Position',[0 0 1 1],'Color','none','XColor','none','YColor','none');
        
        
        if length(size(img))>2        
            text(axtext,.15,rowslice,'slice', 'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom')
%             text(axtext,.92,rowslice,num2str(s), 'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom')
        end
        
        % label sliders
        text(axtext,.15,rowthresh,'threshold', 'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom')
        text(axtext,.15,rowsize,'size filt', 'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom')
        text(axtext,.15,rowsmooth,'size smooth', 'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom')
        text(axtext,.45,rowfancythresh,'dilate size', 'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom')
        text(axtext,.15,rowsat,'saturation', 'Units','normalized','HorizontalAlignment','right','VerticalAlignment','bottom')
        
        % report slider value
%         if numel(T)>1
%             text(axtext,.92,rowthresh,'adapt', 'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom')
%         else
%             text(axtext,.92,rowthresh,num2str(T*maxI*sat), 'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom')
%         end
        text(axtext,.92,rowsize,num2str(sizefilt), 'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom')
        text(axtext,.65,rowsmooth,num2str(sizesmooth), 'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom')
        text(axtext,.7,rowfancythresh,num2str(dilateadapt), 'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom')
    end
end

