%%
figure()
plot(pax_fwhm.idx.clean_upperboundary)

%%
try
roilist.addroi(roianalysis.preprocessed_ch1,'svdtest','rectangle');
catch ME
    disp(ME.message)
    roilist.modifyroi(roianalysis.preprocessed_ch1,'svdtest');
end
%%
roilist.applyvertices(primary_datastruct.stackch1,'svdtest');
testarr = ans;
%%
figure()
imagesc(testarr(:,:,1));

testarr = testarr(:, :, 1:2000);

%%
signal = reshape(testarr, [size(testarr, 1) * size(testarr, 2), size(testarr, 3)]);
[u, s, v] = svd(gpuArray(signal), "econ");

%%
corrMat = corrcoef(u);
figure;
imagesc(corrMat);
colormap jet;
% clim([0 0.7]);
colorbar;