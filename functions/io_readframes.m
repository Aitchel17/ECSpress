function [imgs] = io_readframes(mobj,imgch,xrange,yrange,zrange)
lenx = xrange(2)-xrange(1)+1;
leny = yrange(2)-yrange(1)+1;
lenz = zrange(2)-zrange(1)+1;

imgs = zeros(leny,lenx,lenz);

h = waitbar(0,'Buffering Frames...');

for z = zrange(1):zrange(2)
	waitbar((z-zrange(1))/lenz, h, ...
		sprintf(['Loading:' '%d /' num2str(zrange(2)-zrange(1)+1)], z-zrange(1)));
    img = mobj.ReadFrame(imgch,z)';
	imgs(:,:,z-zrange(1)+1)=img(yrange(1):yrange(2),xrange(1):xrange(2)); % width come first and height come next,  
end


delete(h)
end

