function [trackTable, finalPts] = klt_track_reseed(imgstack, opt)
%KLT_TRACK_RESEED  KLT point tracking with drift reset and re-seeding.
%   Drift is judged by an anchored NCC: a track's current patch is compared with
%   the patch it had when it was born (corr2 against its seed), never with the
%   previous frame. A track that breaks is dropped and a new point is seeded
%   elsewhere; survivors keep their integer id across resets.
%
% IN   imgstack   H x W x N grayscale-able stack
%      opt        name-value tracker settings, defaults in the arguments block;
%                 opt.nccThr is the drift threshold (default 0.5)
% OUT  trackTable long-format table [frame id x y] (px), one row per tracked point
%      finalPts   struct with .id and .xy of the last active set

arguments
    imgstack
    opt.targetN         double = 50     % target number of active tracks
    opt.patchR          double = 15     % anchor patch radius (px)
    opt.nccThr          double = 0.5     % anchored-NCC drift threshold (main knob)
    opt.scoreThr        double = 0.4     % vision.PointTracker score threshold
    opt.maxDisp         double = 12      % max displacement per frame (px)
    opt.minReseedDist   double = 15      % keep re-seeds this far from survivors (px)
    opt.minEigenQuality double = 0.05    % detectMinEigenFeatures MinQuality
    opt.maxBidiErr      double = 2       % vision.PointTracker MaxBidirectionalError (px)
    opt.blockSize       double = [31 31]
    opt.pyrLevels       double = 3
end

% 0. Setup
N = size(imgstack, 3);
frameCells = cell(N, 1);

% 1. Seed on frame 1 and initialize the tracker
% 1.1 Grayscale the frame
I = mat2gray(imgstack(:, :, 1));
% 1.2 Detect strong, border-safe features, capped at targetN
xy = local_detect(I, opt.targetN, zeros(0, 2), opt);
% 1.3 Create and initialize the tracker
tracker = vision.PointTracker('MaxBidirectionalError', opt.maxBidiErr, ...
    'NumPyramidLevels', opt.pyrLevels, 'BlockSize', opt.blockSize);
initialize(tracker, xy, I);
% 1.4 Active-track bookkeeping, kept in the tracker's point order
nextId = 1;
act = local_makeact(I, xy, opt.patchR, nextId);
nextId = nextId + size(xy, 1);
% 1.5 Record frame 1
frameCells{1} = [ones(size(xy, 1), 1), act.id, xy];
needReinit = false;

% 2. Frame loop
for t = 2:N
    I = mat2gray(imgstack(:, :, t));

    % 2.1 Seed fresh and re-initialize when the tracker holds no valid state
    if needReinit || isempty(act.id)
        xy = local_detect(I, opt.targetN, zeros(0, 2), opt);
        if isempty(xy), frameCells{t} = zeros(0, 4); continue; end
        release(tracker); initialize(tracker, xy, I);
        act = local_makeact(I, xy, opt.patchR, nextId);
        nextId = nextId + size(xy, 1);
        needReinit = false;
        frameCells{t} = [t * ones(size(xy, 1), 1), act.id, xy];
        continue;
    end

    % 2.2 Track step
    [xy, isFound, scores] = tracker(I);

    % 2.3 Break detection
    ncc  = local_ncc(I, xy, act.anchor, opt.patchR);      % current vs seed patch
    disp = hypot(xy(:,1) - act.lastxy(:,1), xy(:,2) - act.lastxy(:,2));
    broke = ~isFound ...                    % lost
          | scores < opt.scoreThr ...       % weak match
          | isnan(ncc) | ncc < opt.nccThr ...% drifted off the seed morphology
          | disp > opt.maxDisp;             % teleport
    keep = ~broke;

    % 2.4 Keep survivors with their id, anchor and age
    xy_keep = xy(keep, :);
    act = local_subset(act, keep);
    act.age = act.age + 1;

    % 2.5 Re-seed to refill targetN, clear of the survivors
    nNeed = opt.targetN - size(xy_keep, 1);
    if nNeed > 0
        new_xy = local_detect(I, nNeed, xy_keep, opt);
        if ~isempty(new_xy)
            newAct  = local_makeact(I, new_xy, opt.patchR, nextId);
            nextId  = nextId + size(new_xy, 1);
            act     = local_catact(act, newAct);
            xy_keep = [xy_keep; new_xy]; %#ok<AGROW>
        end
    end

    % 2.6 Push the combined set to the tracker; setPoints makes frame I the new
    %     template source, and an empty set defers to a re-init next frame
    if isempty(xy_keep)
        release(tracker); needReinit = true;
        act = local_emptyact();
        frameCells{t} = zeros(0, 4);
        continue;
    end
    setPoints(tracker, xy_keep);
    act.lastxy = xy_keep;

    % 2.7 Record the frame
    frameCells{t} = [t * ones(size(xy_keep, 1), 1), act.id, xy_keep];
end

% 3. Pack the output
allrows = vertcat(frameCells{:});
trackTable = array2table(allrows, 'VariableNames', {'frame', 'id', 'x', 'y'});
finalPts = struct('id', act.id, 'xy', act.lastxy);

end


% ---------------------------------------------------------------------------

function xy = local_detect(I, nWanted, excludeXY, opt)
%LOCAL_DETECT  Strongest border-safe features away from excludeXY.
%   Returns at most nWanted points as [x y], each at least opt.minReseedDist px
%   from every row of excludeXY.
xy = zeros(0, 2);
if nWanted <= 0, return; end
pts = detectMinEigenFeatures(I, 'MinQuality', opt.minEigenQuality, 'FilterSize', 3);
loc = pts.Location;  met = pts.Metric;
R = opt.patchR;  [H, W] = size(I);
inb = loc(:,1) > R & loc(:,1) <= W-R & loc(:,2) > R & loc(:,2) <= H-R;  % full-patch safe
loc = loc(inb, :);  met = met(inb);
if ~isempty(excludeXY) && ~isempty(loc)
    d = pdist2(loc, excludeXY);
    far = all(d > opt.minReseedDist, 2);
    loc = loc(far, :);  met = met(far);
end
[~, ord] = sort(met, 'descend');
loc = loc(ord, :);
xy = loc(1:min(nWanted, size(loc, 1)), :);
end

function act = local_makeact(I, xy, R, startId)
%LOCAL_MAKEACT  Active-track struct for points xy, ids counting up from startId.
K = size(xy, 1);
act.id     = (startId:startId+K-1).';
act.anchor = local_patches(I, xy, R);
act.age    = ones(K, 1);
act.lastxy = xy;
end

function act = local_emptyact()
%LOCAL_EMPTYACT  Active-track struct with no tracks.
act.id = zeros(0,1);  act.anchor = {};  act.age = zeros(0,1);  act.lastxy = zeros(0,2);
end

function P = local_patches(I, xy, R)
%LOCAL_PATCHES  Fixed-size (2R+1) square seed patch per point.
%   xy must be border-safe; local_detect guarantees that.
K = size(xy, 1);  P = cell(K, 1);
for k = 1:K
    c = round(xy(k, :));
    P{k} = I(c(2)-R:c(2)+R, c(1)-R:c(1)+R);
end
end

function ncc = local_ncc(I, xy, anchors, R)
%LOCAL_NCC  Correlate each current patch with its own seed patch.
%   NaN where the patch would run off the image.
[H, W] = size(I);  K = size(xy, 1);  ncc = nan(K, 1);
for k = 1:K
    c = round(xy(k, :));
    if c(1) > R && c(1) <= W-R && c(2) > R && c(2) <= H-R
        patch = I(c(2)-R:c(2)+R, c(1)-R:c(1)+R);
        ncc(k) = corr2(patch, anchors{k});
    end
end
end

function act = local_subset(act, keep)
%LOCAL_SUBSET  Keep the tracks selected by the logical mask keep.
act.id     = act.id(keep);
act.anchor = act.anchor(keep);
act.age    = act.age(keep);
act.lastxy = act.lastxy(keep, :);
end

function a = local_catact(a, b)
%LOCAL_CATACT  Append active-track struct b after a.
a.id     = [a.id;     b.id];
a.anchor = [a.anchor; b.anchor];
a.age    = [a.age;    b.age];
a.lastxy = [a.lastxy; b.lastxy];
end
