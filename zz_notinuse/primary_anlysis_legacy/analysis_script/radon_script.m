
%%
roilist = roi_handle(fullfile(directories.save_dir,"roilist.mat"));
%%
roilist.addroi(primary_datastruct.stackch2,'radon','rectangle')
%%

mask = roilist.getmask('radon');
%%
holdstack = roilist.applyvertices(primary_datastruct.stackch1,'radon');

%%
sliceViewer(holdstack)
%%
radon_result = analyze_radon(holdstack);

%%
figure()
imagesc(squeeze(radon_result.tirs(:,:,180))==1|squeeze(radon_result.tirs(:,:,180))==3)

%%
figure()
imagesc(squeeze(radon_result.tirs(:,:,1314))==1|squeeze(radon_result.tirs(:,:,1314))==3)

%%

irtdmask = radon_result.irtd>0.03;
%%
sliceViewer(irtdmask)

%%
resize_holdstack = zeros(size(irtdmask));
%%
rgb_stack(diff_coord(1):end-diff_coord(1)-1,diff_coord(2):end-diff_coord(2)-1,:) = holdstack;
%%
x = rgb_stack(diff_coord(1):end,diff_coord(2):end,:,2);
%%
util_checkstack(irtdmask)
