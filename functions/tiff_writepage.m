function tiff_writepage(handle, tags, page, frame)
%TIFF_WRITEPAGE  one page onto an open TIFF
%   Caller: io_postsavetiff
%   A copy of writepage in 00_mdfExtractor's mdf.m; keep the two the same
%
% IN   handle  1 x 1 Tiff      open for writing
%      tags    1 x 1 struct    set on every page, size included
%      page    1 x 1 double    which page of the FILE this is, 1-based
%      frame   H x W uint16
%
%   The directory goes in FRONT of every page but the first, so nothing needs to
%   know how many pages are coming and the file ends without a trailing empty one.
    if page > 1
        handle.writeDirectory();
    end
    handle.setTag(tags);
    handle.write(frame);
end
