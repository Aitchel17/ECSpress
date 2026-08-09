function [boundary_idx,clean_kymograph] = clean_thresholdedkymograph(thresholded, upmod, debug_ax)
%%
x_length = size(thresholded,2);
clean_kymograph = imfill(thresholded,'holes');
boundary_idx = zeros([1,x_length]);

nhood = zeros([7,3]);
if upmod
    nhood(1:2,2) =1; % up
else
    nhood(end-1:end,2) =1; % down
end
clean_kymograph = imerode(clean_kymograph,nhood);
clean_kymograph = imdilate(clean_kymograph,nhood);
clean_kymograph(1:2,:) = 0; % remove up dilation
clean_kymograph(end-1:end,:) = 0; % remove bottom dilation

%%
zero_column = sum(clean_kymograph,1);
zero_column = zero_column ==0;
%%
clean_kymograph(:,zero_column) = thresholded(:,zero_column);
%%
zeros_array = sum(clean_kymograph,1);
%%

%%
%%
cla(debug_ax);
%plot(thresholded(:,8))
% %%
imagesc(debug_ax,clean_kymograph)
%%
for c_idx = 1: x_length
    szfilt.x = clean_kymograph(:,c_idx);
    szfilt.dx =diff([0;szfilt.x;0]);
    szfilt.sx = find(szfilt.dx == 1);
    szfilt.ex = find(szfilt.dx == -1);
    szfilt.lx = szfilt.ex -szfilt.sx;

    if isempty(szfilt.sx)     %no half-max crossing  Mark as NaN (missing) instead
        boundary_idx(c_idx) = NaN;
        continue;
    end

    szfilt.y = ~szfilt.x;
    szfilt.dy =diff([0;szfilt.y;0]);
    szfilt.sy = find(szfilt.dy == 1);
    szfilt.ey = find(szfilt.dy == -1);
    szfilt.ly = szfilt.ey -szfilt.sy;
    [~,szfilt.lxmaxszidx] = max(szfilt.lx);
    %%
    if upmod
        szfilt.count = 1;% up
        szfilt.ly = szfilt.ly(2:end); % up
        while szfilt.count < szfilt.lxmaxszidx
            if szfilt.lx(szfilt.count) < szfilt.ly(szfilt.count)
                % disp(szfilt.lx(szfilt.count))
                szfilt.sx(1:szfilt.count) = Inf;
                szfilt.ex(1:szfilt.count) = Inf;
            end
            szfilt.count=szfilt.count+1;
        end
        boundary_idx(c_idx) = min(szfilt.sx);
    else
        szfilt.count = numel(szfilt.lx); % down
        szfilt.ly = szfilt.ly(1:end-1); % down
        while szfilt.count > szfilt.lxmaxszidx
            if szfilt.lx(szfilt.count) < szfilt.ly(szfilt.count)
                %disp(szfilt.lx(szfilt.count))
                szfilt.sx(szfilt.count:end) = 0;
                szfilt.ex(szfilt.count:end) = 0;
            end
            szfilt.count=szfilt.count-1;
        end
        boundary_idx(c_idx) = max(szfilt.ex);
    end

end
%%
%%
cla(debug_ax);
%%
imagesc(debug_ax,thresholded)
%%
imagesc(debug_ax,clean_kymograph)
%%
end

