savepath = 'E:\OneDrive - The Pennsylvania State University\2023_ECSpress\01_primary_analysis\01_pvscollapse';
vivo_50um1 = mdfExtractLoader();
vivo_50um1.analog = vivo_50um1.loadanalog;
vivo_50um1.stackch1 = vivo_50um1.loadstack("ch1");
vivo_50um1.stackch2 = vivo_50um1.loadstack("ch2");
%%
tmp.gpstack =  analyze_grouproject(vivo_50um1.stackch2,3,"mean");
tmp.medcsf = medfilt3(tmp.gpstack,[1 1 5]);
gausscsf= imgaussfilt(tmp.medcsf,1);
tmp.gpstack = analyze_grouproject(vivo_50um1.stackch1,3,"mean");
tmp.medbv = medfilt3(tmp.gpstack,[1 1 5]);
gaussbv= imgaussfilt(tmp.medbv,1);

%%
roi1 = roi(gaussbv,'pax','line',vivo_50um1);
%%
roi1 = roi1.modifyroi(gausscsf,'pax');
%%
roi1 = roi1.addroi(gaussbv,'sax');
%%
roi_rectangle_polygon(gaussbv,'rectangle');
%%

mdf_rectangle_polygon(gaussbv,'rectangle')
%%
sliceViewer(vivo_50um1.stackch1)


%% primary axis calculation
pax = struct();
pax.rc_bv=analyze_affine_rotate(gaussbv,roi1.vertices.pax(1:2,:), roi1.vertices.pax(3,1));
pax.rc_csf=analyze_affine_rotate(gausscsf,roi1.vertices.pax(1:2,:), roi1.vertices.pax(3,1));
pax.csf_kymograph = squeeze(sum(pax.rc_csf,1));
pax.bv_kymograph = squeeze(sum(pax.rc_bv,1));
[pax.bv_idx, pax.bv_kymomask] = analyze_fwhm(pax.bv_kymograph,0.5,10);
pax.bv_fwhm = pax.bv_idx.lowerboundary_idx-pax.bv_idx.upperboundary_idx;
% reconstruction
tmp.v_thr = repmat(pax.bv_kymomask.upline+pax.bv_kymomask.downline,[1,1,size(pax.rc_bv,1)]);
tmp.v_thr = permute(tmp.v_thr,[3,1,2]);
pax.bvborder = analyze_affine_reverse(tmp.v_thr,size(gaussbv),roi1.vertices.pax(1:2,:));
% normalize projection
tmp.kgph = pax.bv_kymograph;
tmp.kgph = tmp.kgph-min(tmp.kgph,[],1);
pax.norm_kymograph = tmp.kgph./max(tmp.kgph,[], 1);
%%
util_checkstack(~pax.bvborder.*gaussbv) % check reconstructed
%%
figure(Name='normalized vessel kymograph')
imagesc(pax.norm_kymograph)


%% csf processing
[pax.csf_idx, pax.csf_kymomask] = analyze_csfoutter(pax.csf_kymograph,pax.bv_idx.upperboundary_idx,pax.bv_idx.lowerboundary_idx,0.5);
tmp.v_thr = repmat(pax.csf_kymomask.boundary,[1,1,size(pax.rc_bv,1)]);
tmp.v_thr = permute(tmp.v_thr,[3,1,2]);
pax.csfborder = analyze_affine_reverse(tmp.v_thr,size(gaussbv),roi1.vertices.pax(1:2,:));
%%
util_checkstack((~pax.csfborder+~pax.bvborder).*gausscsf)
%%
util_checkstack(~csf_reconstructed.*gausscsf)
%% Secondary axis calculation

roi1 = roi1.modifyroi(gausscsf,'sax');
%%
sax = struct();
sax.rc_bv=analyze_affine_rotate(gaussbv,roi1.vertices.sax(1:2,:), roi1.vertices.sax(3,1));
sax.rc_csf=analyze_affine_rotate(gausscsf,roi1.vertices.sax(1:2,:), roi1.vertices.sax(3,1));
sax.csf_kymograph = squeeze(sum(sax.rc_csf,1));
sax.bv_kymograph = squeeze(sum(sax.rc_bv,1));
[sax.bv_idx, sax.bv_kymomask] = analyze_fwhm(sax.bv_kymograph,0.5,10);
sax.bv_fwhm = sax.bv_idx.lowerboundary_idx-sax.bv_idx.upperboundary_idx;
% reconstruction
tmp.v_thr = repmat(sax.bv_kymomask.upline+sax.bv_kymomask.downline,[1,1,size(sax.rc_bv,1)]);
tmp.v_thr = permute(tmp.v_thr,[3,1,2]);
sax.bvborder = analyze_affine_reverse(tmp.v_thr,size(gaussbv),roi1.vertices.sax(1:2,:));
% normalize projection
tmp.kgph = sax.bv_kymograph;
tmp.kgph = tmp.kgph-min(tmp.kgph,[],1);
sax.norm_kymograph = tmp.kgph./max(tmp.kgph,[], 1);
%%
[sax.csf_idx, sax.csf_kymomask] = analyze_csfoutter(sax.csf_kymograph,sax.bv_idx.upperboundary_idx,sax.bv_idx.lowerboundary_idx,0.5);
tmp.v_thr = repmat(sax.csf_kymomask.boundary,[1,1,size(sax.rc_bv,1)]);
tmp.v_thr = permute(tmp.v_thr,[3,1,2]);
sax.csfborder = analyze_affine_reverse(tmp.v_thr,size(gaussbv),roi1.vertices.sax(1:2,:));
%%
figure()
plot(medfilt1(pax.csf_idx.upboundary,5))
hold on
plot(medfilt1(pax.bv_idx.upperboundary_idx,5))
%%
figure()
plot(medfilt1(pax.csf_idx.downboundary,5))
hold on
plot(medfilt1(pax.bv_idx.lowerboundary_idx,5))
%%
figure()
plot(medfilt1(sax.csf_idx.upboundary,5))
hold on
plot(medfilt1(sax.bv_idx.upperboundary_idx,5))
%%
figure()
plot(medfilt1(sax.csf_idx.downboundary,5))
hold on
plot(medfilt1(sax.bv_idx.lowerboundary_idx,5))
%%
util_checkstack(~sax.csfborder.*gaussbv) % check reconstructed
util_checkstack((~sax.csfborder+~sax.bvborder+~pax.csfborder+~pax.bvborder).*gausscsf) % check reconstructed


%%
figure(Name='normalized vessel kymograph')
imagesc(sax.norm_kymograph)
