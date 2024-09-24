[info, tpsm_info, mode_info, mobj] = io_info();

%%
%generate chunk
imgch =1;

startIdx   = 1;
endIdx     = tpsm_info.fcount;
imgs = zeros(tpsm_info.fh,tpsm_info.fw,tpsm_info.fcount);

h = waitbar(0,'Buffering Frames...');

for idx= startIdx:endIdx
	waitbar((idx-startIdx)/tpsm_info.fcount, h, ...
		sprintf(['Loading:' '%d /' num2str(endIdx-startIdx+1)], idx-startIdx));
    img = mobj.ReadFrame(imgch,1+idx-startIdx)';
	imgs(:,:,idx)=img; % width come first and height come next,  
end


delete(h)

%%
imshow(imgs(:,:,330)',[0,400])

%%
analogsample = str2double(mobj.ReadParameter('Analog Sample Count'));

%%
analog = mobj.ReadAnalog(5,analogsample,0);
%%

nanarray = imgs == -2048;


%%
plot(analog)
%
for i = 0:1:7
    strcat('mode_info.analog'+string(i))
end
%%
imshow(imgs(:,:,5)');