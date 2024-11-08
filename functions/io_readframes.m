function zstack = io_readframes(mobj,imgch,zrange)

lenz = zrange(2)-zrange(1)+1;
tic
y = arrayfun(@(z) loadimg(z), (zrange(1):zrange(2)), 'un',0);
zstack = cat(3,y{:});
toc
    function frame = loadimg(z)
        fprintf("\r Reading frames from .mdf %d%%", round(z*100/lenz));
        frame = mobj.ReadFrame(imgch,z)';
    end
end

