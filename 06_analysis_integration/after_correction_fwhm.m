%% analysis window initialization
window1=figure(Name='window1',NumberTitle='off');
window2=figure(Name='window2',NumberTitle='off');
% window3=figure(Name='window3',NumberTitle='off');
% window4=figure(Name='window4',NumberTitle='off');

%% gaussian fit full width half max
test.testedlist = {'HQL072_250325_003',linedata_struct.HQL072_250325_003.infomap('Comments');
  'HQL072_250331_005',  linedata_struct.HQL072_250331_005.infomap('Comments');
  'HQL072_250331_003', linedata_struct.HQL072_250331_003.infomap('Comments');
  'HQL072_250331_008', linedata_struct.HQL072_250331_008.infomap('Comments'); };
%%
test.testedlist('HQL080_250722_002')
%%
clee = color_lee;
test.snippet_fname = 'HQL080_whiskerb_250722_002';
disp(linedata_struct.(test.snippet_fname).infomap('Comments'))

test.pax = linedata_struct.(test.snippet_fname).fwhmline.pax;
test.kgph_bv = test.pax.kymograph.kgph_bv;
test.kgph_csf = test.pax.kymograph.kgph_csf;
%1 using ceterbvidx, seperate bottom half and upperhalf
%2 using minimum position identification might screwed <<< lets try
%without min position 
test.centerbvidx = ceil((test.pax.idx.bv_lowerboundary+test.pax.idx.bv_upperboundary)/2);
test.rowidx = test.pax.mask.rowidx;

tmp.upbv = test.pax.idx.bv_upperboundary;
tmp.upcell = test.pax.idx.pvs_upshadeloc;
tmp.uppvs = test.pax.idx.pvs_upperboundary;
%
tmp.upbv = medfilt1(tmp.upbv,20);
tmp.upcell = medfilt1(tmp.upcell,20);
tmp.uppvs = medfilt1(tmp.uppvs,20);
%
tmp.downbv = test.pax.idx.bv_lowerboundary;
tmp.downcell = test.pax.idx.pvs_downshadeloc;
tmp.downpvs = test.pax.idx.pvs_lowerboundary;
%
tmp.downbv = medfilt1(tmp.downbv,20);
tmp.downcell = medfilt1(tmp.downcell,20);
tmp.downpvs = medfilt1(tmp.downpvs,20);
%
tmp.downbv = tmp.downbv(2:end);
tmp.downcell = tmp.downcell(2:end);
tmp.downpvs = tmp.downpvs(2:end);

tmp.upbv = tmp.upbv(2:end);
tmp.upcell = tmp.upcell(2:end);
tmp.uppvs = tmp.uppvs(2:end);
%
test.bv_thickness = tmp.downbv - tmp.upbv;
test.uppvs_thickness = tmp.upbv - tmp.uppvs;
test.downpvs_thickness = tmp.downpvs - tmp.downbv;
test.uppvs_changes = test.uppvs_thickness-min(test.uppvs_thickness);
test.downpvs_changes = test.downpvs_thickness-min(test.downpvs_thickness);
test.bv_changes = test.bv_thickness-min(test.bv_thickness);


%%
figure(window1)
hold off
plot(tmp.upbv)
hold on
%plot(tmp.upcell)
plot(tmp.uppvs)
plot(tmp.upbv-tmp.uppvs)

figure(window2)
hold off
plot(tmp.downbv)
hold on
%plot(tmp.downcell)
plot(tmp.downpvs)
plot(tmp.downpvs-tmp.downbv)

%%



%%
figure()
plot(test.uppvs_changes)
hold on
plot(test.downpvs_changes)
plot(test.bv_changes)
%%

%%
plot_2dhistdistribution(test.uppvs_changes, test.bv_changes,1/0.57,'down pvs',clee.gray_g)
%%

h2frame = getframe(gcf);
h2img = frame2im(h2frame);
%%
figure()
imshow(h2img)
%% first data analysis
crop = [];

if isempty(crop)
    x_data = test.uppvs_changes;
    y_data = test.bv_changes;
else
    x_data = test.uppvs_changes(crop(1):crop(2));
    y_data = test.bv_changes(crop(1):crop(2));
end

scale = 1/0.29;
figtitle = 'down pvs';
cmap = clee.gray_g;
% Create 2D histogram with fixed bin width
scaled_x = x_data/scale;
scaled_y = y_data/scale;
[x_counts, xedges, yedges] = histcounts2(x_data/scale, y_data/scale, 'BinWidth', [1/scale, 1/scale],'Normalization','percentage');
% Get bin centers
x_centers = (xedges(1:end-1) + xedges(2:end)) / 2;
bv_centers = (yedges(1:end-1) + yedges(2:end)) / 2;

% Apply log scale
log_values = log10(x_counts + 1); % Avoid log(0)

% second data analysis
if isempty(crop)
    x_data = test.downpvs_changes;
    y_data = test.bv_changes;
else
    x_data = test.downpvs_changes(crop(1):crop(2));
    y_data = test.bv_changes(crop(1):crop(2));
end

% scale = 1/0.57;
figtitle = 'down pvs';
cmap = clee.gray_g;
% Create 2D histogram with fixed bin width
scaled_x = x_data/scale;
scaled_y = y_data/scale;
[x_counts, xedges, yedges] = histcounts2(x_data/scale, y_data/scale, 'BinWidth', [1/scale, 1/scale],'Normalization','percentage');
% Get bin centers
x_centers2 = (xedges(1:end-1) + xedges(2:end)) / 2;
bv_centers = (yedges(1:end-1) + yedges(2:end)) / 2;

% Apply log scale
log_values2 = log10(x_counts + 1); % Avoid log(0)

% third data analysis
if isempty(crop)
    x_data = test.downpvs_changes + test.uppvs_changes;
    y_data = test.bv_changes;
else
    x_data = test.downpvs_changes(crop(1):crop(2)) + test.uppvs_changes(crop(1):crop(2));
    y_data = test.bv_changes(crop(1):crop(2));
end
x_data = x_data-min(x_data);

% scale = 1/0.57;
figtitle = 'down pvs';
cmap = clee.gray_g;
% Create 2D histogram with fixed bin width
scaled_x = x_data/scale;
scaled_y = y_data/scale;
[x_counts, xedges, yedges] = histcounts2(x_data/scale, y_data/scale, 'BinWidth', [1/scale, 1/scale],'Normalization','percentage');
% Get bin centers
x_centers3 = (xedges(1:end-1) + xedges(2:end)) / 2;
bv_centers = (yedges(1:end-1) + yedges(2:end)) / 2;

% Apply log scale
log_values3 = log10(x_counts + 1); % Avoid log(0)

%
x_centerarray = {x_centers3,x_centers2,x_centers};
[max_pvschange,maxloc] = max([length(x_centerarray{1}),length(x_centerarray{2}),length(x_centerarray{3})]);
max_centerarray = x_centerarray{maxloc};
%%
spattern = "flat";
resize_logvalue = zeros([max_pvschange,length(bv_centers)]);
resize_logvalue(1:size(log_values3,1),1:size(log_values3,2))= log_values3;

figure()
s = pcolor(bv_centers, max_centerarray, resize_logvalue);
s.FaceColor = spattern;
set(s,'EdgeColor','none');
set(gca, 'YDir', 'normal')
axis square
axis off
colormap(clee.gray_g)
h3frame = getframe(gca);
h3img = frame2im(h3frame);



resize_logvalue = zeros([max_pvschange,length(bv_centers)]);
resize_logvalue(1:size(log_values2,1),1:size(log_values2,2))= log_values2;
figure()
s = pcolor(bv_centers, max_centerarray, resize_logvalue);
s.FaceColor = spattern;
set(s,'EdgeColor','none');
set(gca, 'YDir', 'normal')
axis square
axis off


colormap(clee.red_g)
%
resize_logvalue2 = zeros(size(log_values));
resize_logvalue2(1:size(log_values2,1),1:size(log_values2,2)) = log_values2;
h2frame = getframe(gca);
h2img = frame2im(h2frame);


resize_logvalue = zeros([max_pvschange,length(bv_centers)]);
resize_logvalue(1:size(log_values,1),1:size(log_values,2))= log_values;
figure()
s = pcolor(bv_centers, max_centerarray, resize_logvalue);
axis square;
axis off

colormap(clee.green_g)
s.AlphaData = 0.5;
s.FaceColor = spattern;
set(s,'EdgeColor','none');
h1frame = getframe(gca);
h1img = frame2im(h1frame);

figure()
imshow(h1img+h2img+h3img)
%%
figure()
imshow(h1frame.cdata)

%%
ax2.Visible = 'off'; 
ax2.XTick = []; 
ax2.YTick = []; 
colormap(ax1,clee.red_g) 
colormap(ax2,clee.gray_g) 
%%
xlabel('bv thickness / mean')
ylabel('downpvs thickness')
title(figtitle)

















%%
%%
figure()
scatter(test.bv_thickness,test.uppvs_thickness)

%%
imagesc(histogram2)

%%
y_cell = {};
%%
y_cell = [y_cell,1];
%%


%%
figure()
plot(test.bv(2:end))
%%
figure()
plot(test.downpvs(2:end))

%%
updown = 'uppvs_thickness';
x = test.bv(2:end);
y = test.(updown)(2:end);
x_unique = unique(x);
mean_y = zeros(size(x_unique));
ci_y = zeros(size(x_unique));  % 95% CI = 1.96*SEM
y_cell = {};
for i = 1:length(x_unique)
    xi = x_unique(i);
    yi = y(x == xi);
    y_cell = [y_cell,yi];
    mean_y(i) = mean(yi);
    SEM = std(yi) / sqrt(length(yi)); 
    tval = tinv(0.975, length(yi) - 1);
    ci_y(i) = tval * SEM;              
end

figure
errorbar(x_unique, mean_y, ci_y, 'o-', 'Color', 'k', ...
    'MarkerFaceColor', 'r', 'LineWidth', 1.5)
hold on
xlabel('BV diameter')
ylabel('up PVS thickness')
title(['Mean ± 95% CI of ', updown ,' for each bv'])
%%




%% box plot
figure(window1)
hold off
boxplot(test.uppvs(2:end),test.bv(2:end))
%
figure(window2)
hold off
boxplot(test.downpvs(2:end),test.bv(2:end))


%%
test.kgph_csf_up = test.kgph_csf;
test.kgph_csf_up(test.rowidx<=test.centerbvidx) = NaN;

%%
test.kgph_csf_down = test.kgph_csf;
test.kgph_csf_down(test.rowidx>=test.centerbvidx) = NaN;
%%
figure(window1)
imagesc(test.kgph_csf_up)

figure(window2)
imagesc(test.kgph_csf_down)

figure(window3)
imagesc(test.kgph_bv)

figure(window4)
plot(test.kgph_csf_down(:,100:150))
%%
tmp.norm_kgph_csf = test.kgph_csf./max(test.kgph_csf,[],1);

%%
figure(window4)
for i = 1:1000
    oneline = test.kgph_csf_up(:,i);
    x = (1:length(oneline))';
    
    % Create logical exclusion mask
    excludeMask = isnan(oneline);
    
    % Optional: replace NaNs with 0 or any filler just to visualize, not needed for fitting
    
    % Fit Gaussian (2-peak) model, excluding NaNs
    fit_result = fit(x, oneline, 'gauss3', 'Exclude', excludeMask);
    %
    hold off
    plot(oneline)
    plot(x, oneline, 'k-', 'LineWidth', 1.5); hold on;
    
    % Overlay fit result
    plot(fit_result, 'r-');  % this automatically uses same x range
    legend('Original Data', 'Gaussian Fit');
    xlabel('Y (Distance)');
    ylabel('Intensity');
    title('Gaussian Fit Overlayed on Original Line');
    grid on;
end
%%
%% Load data
kg = test.kgph_csf_down; % [distance x time]
y_axis = (1:size(kg,1))';

% Prepare Gaussian fit parameters
fit_type = 'gauss2'; % use single gaussian (use 'gauss2' if preferred)
kg_fit = zeros(size(kg)); % Store Gaussian-fitted kymograph

% Loop over each column, fit Gaussian
for col = 1:100
    disp(col)
    intensity = kg(:, col);
    
    excludeMask = isnan(intensity);
    
    % Gaussian fitting with exclusion
    fit_result = fit(y_axis, intensity, fit_type, 'Exclude', excludeMask);
        
        % Store parameters (A, mu, sigma)
    coeffs = coeffvalues(fit_result);
        
        % Reconstruct fitted Gaussian curve
    fitted_curve = feval(fit_result, y_axis);
    kg_fit(:,col) = fitted_curve;

end

%%

figure(window4)
hold off
plot(kg(:,60))
hold on
% plot(fit_result)
plot( kg_fit(:,60))
%% Visualization: Original vs. Gaussian-fitted kymograph
figure;

subplot(2,1,1);
imagesc(kg);
title('Original Kymograph');
xlabel('Time');
ylabel('Distance');
colormap(turbo);
colorbar;
axis tight;

subplot(2,1,2);
imagesc(kg_fit);
title('Gaussian-Fitted Kymograph');
xlabel('Time');
ylabel('Distance');
colormap(turbo);
colorbar;
axis tight;
