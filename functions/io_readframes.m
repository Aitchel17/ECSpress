function [imgs] = io_readframes(mobj,imgch,xrange,yrange,zrange)
lenx = xrange(2)-xrange(1)+1;
leny = yrange(2)-yrange(1)+1;
lenz = zrange(2)-zrange(1)+1;

imgs = zeros(leny,lenx,lenz);

h = waitbar(0,'Buffering Frames...');




loadingnumber = 10;
zmod = mod(lenz,loadingnumber);

disp('start load leftover frame')
disp(lenz)
disp(zmod)
if zmod ~= 0
    for z = zrange(1):zrange(1)+zmod-1
	    waitbar((z-zrange(1))/lenz, h, ...
	    sprintf(['Loading:' '%d /' num2str(zrange(2)-zrange(1)+1)], z-zrange(1)));
        img = mobj.ReadFrame(imgch,z)';
        imgs(:,:,z-zrange(1)+1)=img(yrange(1):yrange(2),xrange(1):xrange(2)); % width come first and height come next, 
    end
end
disp(zrange(1))
disp(zrange(1)+zmod-1)
disp('start load 10 frames')
disp(zrange(1)+zmod)
disp(zrange(2)-zmod-loadingnumber+1)
for z = zrange(1)+zmod:loadingnumber:zrange(2)-zmod-loadingnumber+1
	waitbar((z-zrange(1))/lenz, h, ...
	sprintf(['Loading:' '%d /' num2str(zrange(2)-zmod-loadingnumber+1)], z-zrange(1)+zmod+loadingnumber));
    tmpimg = cat(3,mobj.ReadFrame(imgch,z)',mobj.ReadFrame(imgch,z+1)', ...
                mobj.ReadFrame(imgch,z+3)',mobj.ReadFrame(imgch,z+4)', ...
                 mobj.ReadFrame(imgch,z+5)',mobj.ReadFrame(imgch,z+6)', ...
                  mobj.ReadFrame(imgch,z+7)',mobj.ReadFrame(imgch,z+8)', ...
                mobj.ReadFrame(imgch,z+9)',mobj.ReadFrame(imgch,z+10)');
    imgs(:,:,z-zrange(1)+1+zmod:z-zrange(1)+zmod+loadingnumber)=tmpimg(yrange(1):yrange(2),xrange(1):xrange(2),:); % width come first and height come next, 
end
% z - zrange(1) : adjust initial position
%               + zmod : adjust position after leftover
%               

delete(h)
end

