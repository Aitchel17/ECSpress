addpath(genpath(pwd));
sessiondir ='G:\tmp\01_igkltdt\hql086\250909_hql086_sleep\HQL086_sleep250909_005_piv';
session = ECSSession(sessiondir);
session = session.load_primary_results();
sleep_integrate = state_integration(sessiondir);
%todo 2026 04 28
% 1. state transition
% 2. PIV
% 3. Make offset video
% 1. Load data & 3. Load processed data (Integrated via ECSSession)
session.stackch1 = session.loadstack('ch1');
session.stackch2 = session.loadstack('ch2');
% 2. Twophoton data FPS matching & preprocessing
% Note: twophoton_preprocess expects a struct with stackch1/2 and img_param.
% ECSSession object works here as it has these properties.
twophoton_processed = twophoton_preprocess(session);
%%
img_state = state_image(sleep_integrate);
img_state.get_state_indices(twophoton_processed.t_axis,twophoton_processed.outfps);
img_state.param.pixel2um = twophoton_processed.pixel2um;

na_trans = get_stateframes(img_state.state_idx.na_trans,twophoton_processed.ch1);
ra_trans = get_stateframes(img_state.state_idx.ra_trans,twophoton_processed.ch1);
awake = get_stateframes(img_state.state_idx.awake,twophoton_processed.ch1);


of.state = ra_trans{1};
of.pre = piv_preprocess(of.state);

imgStack = of.pre;
I0 = imgStack(:,:,1);
I0 = mat2gray(I0);

points = detectMinEigenFeatures(I0, ...
    'MinQuality', 0.2, ...
    'FilterSize', 3);

xy = points.Location;   % [x y] coordinates

tracker = vision.PointTracker( ...
    'MaxBidirectionalError', 3, ...
    'NumPyramidLevels', 3, ...
    'BlockSize', [31 31]);

tracker.initialize(xy,I0)

% 결과 저장용
nFrames = size(imgStack, 3);
tracks = nan(size(xy,1), 2, nFrames);
valids = false(size(xy,1), nFrames);

tracks(:,:,1) = xy;
valids(:,1) = true;

% 3) 프레임별 추적
for t = 2:nFrames
    It = mat2gray(imgStack(:,:,t));

    [xy_new, isValid] = tracker(It);

    tracks(:,:,t) = xy_new;
    valids(:,t) = isValid;
end

%%

figure;
ax = axes;
frame1 = mat2gray(imgStack(:,:,81));
frame2 = mat2gray(imgStack(:,:,91));
himg = imshow(frame2-frame1);
%clim([-0.2 0.2]);
hold(ax, 'on')
axis(ax, 'image');

% color gradient across time
cmap = turbo(numel(frameRange)); 
% optional: keep plot handles so we can delete/update them
hPts = gobjects(numel(frameRange),1);

for k = 1:numel(frameRange)
    t = frameRange(k);

    % update image
    himg.CData = mat2gray(imgStack(:,:,t));
    %clim([-0.2 0.2]);

    % get valid points
    xy_t = tracks(:,:,t);
    isValid = valids(:,t);

    % plot current frame's points with time-coded color
    hPts(k) = plot(ax, xy_t(isValid,1), xy_t(isValid,2), '.', ...
        'Color', cmap(k,:), ...
        'MarkerSize', 10);

    title(ax, sprintf('KLT points up to frame %d', t));

    drawnow;
    pause(0.7);
end


%%
% frames to show
frameRange = 80:110;

figure;
ax = axes;
hImg = imshow(mat2gray(imgStack(:,:,frameRange(1))), [], 'Parent', ax);
hold(ax, 'on');
axis(ax, 'image');

% color gradient across time
cmap = turbo(numel(frameRange)); 
% optional: keep plot handles so we can delete/update them
hPts = gobjects(numel(frameRange),1);

for k = 1:numel(frameRange)
    t = frameRange(k);

    % update image
    hImg.CData = mat2gray(imgStack(:,:,t));

    % get valid points
    xy_t = tracks(:,:,t);
    isValid = valids(:,t);

    % plot current frame's points with time-coded color
    hPts(k) = plot(ax, xy_t(isValid,1), xy_t(isValid,2), '.', ...
        'Color', cmap(k,:), ...
        'MarkerSize', 10);

    title(ax, sprintf('KLT points up to frame %d', t));

    drawnow;
    if k == 1 
        pause(10);
    end
        pause(0.7);
end


%%

% frames to show
frameRange = 81:111;

figure;
ax = axes;
hImg = imshow(mat2gray(imgStack(:,:,frameRange(1))), [], 'Parent', ax);
hold(ax, 'on');
axis(ax, 'image');

% color gradient across time
cmap = turbo(numel(frameRange));  % or parula, jet, hsv

% optional: keep plot handles so we can delete/update them
hPts = gobjects(numel(frameRange),1);

t =  88;
% update image
hImg.CData = mat2gray(imgStack(:,:,t));

% get valid points
xy_t = tracks(:,:,t);
isValid = valids(:,t);

% plot current frame's points with time-coded color
hPts(k) = plot(ax, xy_t(isValid,1), xy_t(isValid,2), '.', ...
    'Color', cmap(k,:), ...
    'MarkerSize', 10);

title(ax, sprintf('KLT points up to frame %d', t));

drawnow;
pause(0.3);


%%

frameRange = 60:110;
trailLength = 20;  % 최근 10프레임만 표시

figure;
ax = axes;
hImg = imshow(mat2gray(imgStack(:,:,frameRange(1))), [], 'Parent', ax);
hold(ax, 'on');
axis(ax, 'image');

cmap = turbo(trailLength);

for k = 1:numel(frameRange)
    t = frameRange(k);

    hImg.CData = mat2gray(imgStack(:,:,t));

    % delete previous point overlays
    delete(findobj(ax, 'Tag', 'KLTpoints'));

    % show recent trail
    recentIdx = max(1, k-trailLength+1):k;

    for jj = 1:numel(recentIdx)
        kk = recentIdx(jj);
        tt = frameRange(kk);

        xy_tt = tracks(:,:,tt);
        isValid = valids(:,tt);

        plot(ax, xy_tt(isValid,1), xy_tt(isValid,2), '.', ...
            'Color', cmap(jj,:), ...
            'MarkerSize', 10, ...
            'Tag', 'KLTpoints');
    end

    title(ax, sprintf('KLT point trail: frame %d', t));

    drawnow;
    pause(0.8);
end
%%

frameRange = 60:110;
trailLength = 1;  % 최근 10프레임만 표시

figure;
ax = axes;
hImg = imshow(mat2gray(imgStack(:,:,frameRange(1))), [], 'Parent', ax);
hold(ax, 'on');
axis(ax, 'image');

cmap = turbo(trailLength);

t = 60;

    hImg.CData = mat2gray(imgStack(:,:,t));

    % delete previous point overlays
    delete(findobj(ax, 'Tag', 'KLTpoints'));

    % show recent trail
    recentIdx = max(1, k-trailLength+1):k;

    for jj = 1:numel(recentIdx)
        kk = recentIdx(jj);
        tt = frameRange(kk);

        xy_tt = tracks(:,:,tt);
        isValid = valids(:,tt);

        plot(ax, xy_tt(isValid,1), xy_tt(isValid,2), '.', ...
            'Color', cmap(jj,:), ...
            'MarkerSize', 10, ...
            'Tag', 'KLTpoints');

    end
    title(ax, sprintf('KLT point trail: frame %d', t));
