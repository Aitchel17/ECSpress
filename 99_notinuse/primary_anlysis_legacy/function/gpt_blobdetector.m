function [mask, resp, sigmaMap, radii_out] = gpt_blobdetector(I, radii, prc, nScales)
% I: HxW 또는 HxWxT (NaN 무시)
% radii: 비우면([]) 자동 추정, 아니면 [rmin:rstep:rmax] 또는 벡터
% prc: 응답 퍼센타일 임계 (예: 90~99)
% nScales: radii 자동일 때 사용할 스케일 개수(기본 12)
if nargin < 3 || isempty(prc), prc = 95; end
if nargin < 4 || isempty(nScales), nScales = 12; end

is3d = ndims(I)==3; if ~is3d, I = reshape(I,size(I,1),size(I,2),1); end
[H,W,T] = size(I);

% --- 1) radii 자동 결정 (대략적인 범위 추정) ---
if nargin < 2 || isempty(radii)
    % 상위 99%로 씨드 얻고 영역 크기에서 반지름 추정
    It = double(I(:,:,min(T,1)));
    val = It(~isnan(It));
    if isempty(val)
        rMin = 1.2; rMax = min([floor(min(H,W)/3), 30]);    % 완전 기본
    else
        seed = false(H,W); seed(~isnan(It)) = It(~isnan(It)) >= prctile(val, 99);
        CC = bwconncomp(seed, 8);
        if CC.NumObjects >= 3
            A = cellfun(@numel, CC.PixelIdxList);
            rEst = sqrt(A/pi);                    % area→반지름
            rMin = max(1.2, prctile(rEst, 20));   % 과소평가 방지
            rMax = max(rMin*1.5, prctile(rEst, 95)*2); % 약간 여유
        else
            rMin = 1.2; rMax = min([floor(min(H,W)/3), 30]);
        end
    end
    % 로그 간격 스케일
    radii = exp(linspace(log(rMin), log(rMax), nScales));
end
radii_out = radii(:).';
sigmas = radii_out./sqrt(2);   % r ≈ √2·σ

% --- 2) 멀티-스케일 LoG 응답의 최대값/스케일 선택 ---
resp = zeros(H,W,T);
sigmaMap = zeros(H,W,T);
for t = 1:T
    It = double(I(:,:,t));
    valid = ~isnan(It);
    Rmax = -inf(H,W); Smap = zeros(H,W);
    for s = 1:numel(sigmas)
        sg = sigmas(s);
        fsz = max(5, 2*ceil(3*sg)+1);
        h   = fspecial('log', fsz, sg);
        R   = -sg^2 * imfilter(It, h, 'replicate', 'conv');   % 밝은 블롭 +
        R(~valid) = -inf;
        upd = R > Rmax;
        Rmax(upd) = R(upd);
        Smap(upd) = sg;
    end
    resp(:,:,t) = Rmax;
    sigmaMap(:,:,t) = Smap;
end

% --- 3) 퍼센타일 임계로 마스크 생성 (프레임별) ---
mask = false(H,W,T);
for t = 1:T
    r = resp(:,:,t);
    v = r(isfinite(r) & r>0);
    if isempty(v), continue; end
    th = prctile(v, prc);
    mask(:,:,t) = r >= th;
end

if ~is3d
    resp = resp(:,:,1); sigmaMap = sigmaMap(:,:,1); mask = mask(:,:,1);
end
end