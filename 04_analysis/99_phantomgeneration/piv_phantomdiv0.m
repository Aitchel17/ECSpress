%PIV_PHANTOMDIV0  발산이 정확히 0 인 변위장을 실제 프레임에 얹어, PIV 사슬이
%무엇을 되돌려주는지 잰다. Ctrl+Enter 로 셀 단위 실행.
%
%   같은 폴더의 analysis_phantom 과 묻는 것이 다르다.
%     analysis_phantom   이미지를 합성한다. 대칭으로 확장하는 혈관을 그려서 원본을
%                        대신한다 -> FWHM 과 기하 사슬을 시험
%     piv_phantomdiv0    실제 프레임을 알려진 장으로 워프한다 -> PIV 와 발산 사슬
%
%   얹는 장은
%       u = (A/r) * rhat,   A = [(r0 + dD/2)^2 - r0^2] / 2
%   원점 밖에서 발산이 정확히 0 이다. 어느 원을 지나든 같은 면적이 지나가므로
%   (2*pi*r * A/r = 2*pi*A, r 이 지워진다) 임의의 고리에 대한 volume_out 의
%   참값이 0 이고, 파이프라인이 돌려주는 값은 전부 파이프라인 자신의 것이다.
%
%   이득 k(r) 이 반경에 딸리면 발산이 만들어진다. 균일한 이득은 못 만든다.
%   그래서 셀 9 의 k(r) 이 평평한지가 셀 10 의 volume_out 을 결정한다.
%
%   셀 0-1   준비        설정과 세션 확인
%   셀 2-4   정의        기하 -> 진폭 dD -> A
%   셀 5-6   만들기      워프를 얹을 프레임 -> 프레임 쌍
%   셀 7     확인        만든 장을 눈으로 본다. 여기서 안 맞으면 뒤가 무의미
%   셀 8-10  재기        PIV -> 이득 k(r) -> volume_out
%
%   IN   base 워크스페이스의 S / twophoton_processed / event_det / session.
%        piv_integration_testbed 의 셀 0 (session, twophoton_processed),
%        셀 3 (event_det), 셀 5 (S) 가 만든다
%   OUT  워크스페이스에 남는 것. 디스크에 쓰지 않는다
%
%   see CLAUDE_LOG.md

%% 0. 설정
param.condition     = "rem2awake";           % str          이 아래는 조건 하나만
param.n_wedge       = 8;                     % int          등각 분할
param.bin_edges_um  = 0:1.5:40;              % 1 x nB+1 float  벽에서의 거리, um
param.bin_read      = 9;                     % int          벽에서 12.0 - 13.5 um
param.window_sizes  = [40 20; 20 10; 12 6];  % P x 2 int    [창 스텝], 패스마다
param.n_pair        = 5;                     % int          워프를 얹을 프레임 쌍
param.use_gpu       = true;                  % bool
param.min_tri_wedge = 10;                    % int          웨지가 갖춰야 할 삼각형
param.ring_edges_um = 12:4:40;               % 1 x nR+1 float  이득을 재는 고리

dirs.store = '';   % char   claude_swingpiv_run.mat 과 _profile.mat 이 있는 폴더

%% 1. 세션이 올라와 있나
%  깊은 곳에서 실패하는 대신 여기서 말한다
for tmp_v = ["S" "twophoton_processed" "event_det" "session"]
    if ~evalin('base', sprintf('exist(''%s'',''var'')', tmp_v))
        error('piv_phantomdiv0:noSession', ...
            '%s 가 base 에 없다. piv_integration_testbed 의 셀 0, 3, 5 를 먼저.', tmp_v);
    end
end
S                   = evalin('base', 'S');
twophoton_processed = evalin('base', 'twophoton_processed');
event_det           = evalin('base', 'event_det');
session             = evalin('base', 'session');
p2u = twophoton_processed.pixel2um;
[H, W, n_frame] = size(S);
fprintf('프레임 %d x %d x %d,  픽셀 %.4f um\n', H, W, n_frame, p2u);

%% 2. 기하 — 팬텀이 어디를 중심으로 도는가
%  팬텀은 coremask 의 무게중심을 축으로 축대칭이다. bin 은 반경이 아니라 벽에서의
%  거리이므로, 마스크가 원이 아니면 같은 bin 이 방향마다 다른 참 반경에 앉는다
run_store = load(fullfile(dirs.store, 'claude_swingpiv_run.mat'), ...
    'piv_run', 'coremask', 'exclmask');
coremask = run_store.coremask;
exclmask = run_store.exclmask;

vp_geo = vfield_polar(run_store.piv_run(1).xyuv, coremask, p2u, ...
    n_wedge = param.n_wedge, bin_edges_um = param.bin_edges_um, ...
    exclmask = exclmask, gated = true);
centroid_x = vp_geo.centroid(2);
centroid_y = vp_geo.centroid(1);

bv_px = double(session.pax_fwhm.idx.clean_lowerBVboundary(:)).' ...
      - double(session.pax_fwhm.idx.clean_upperBVboundary(:)).';
r0_um = median(bv_px) * p2u / 2;
r_wall_um = sqrt(nnz(coremask)/pi) * p2u;

fprintf('무게중심   x %.1f  y %.1f  px\n', centroid_x, centroid_y);
fprintf('내경 r0    %.3f um     루멘 반지름. 팬텀이 미는 벽\n', r0_um);
fprintf('coremask   %.3f um     등가 반경. bin 은 이 경계에서 잰다\n', r_wall_um);

%% 3. 이 조건의 진폭 dD
%  실제 이벤트들의 직경 변화 중앙값. 팬텀은 이 값 하나에 맞춘다
prof_store = load(fullfile(dirs.store, 'claude_swingpiv_profile.mat'), ...
    'claude_arm', 'claude_all');
row_arm   = string(prof_store.claude_arm);
row_state = string({prof_store.claude_all.state});
row_pol   = string({prof_store.claude_all.pol});
row_dD_px = [prof_store.claude_all.diameter_change];

switch param.condition
    case "rem2awake"
        row_sel = row_state == "rem2awake";
    case "nrem2awake"
        row_sel = row_state == "nrem2awake";
    case "nrem_swing"
        row_sel = row_arm == "swing" & row_state == "nrem" & row_pol == "dilation";
    case "awake_swing"
        row_sel = row_arm == "swing" & row_state == "awake" & row_pol == "dilation";
    otherwise
        error('piv_phantomdiv0:cond', '모르는 조건 %s', param.condition);
end
dD_um = median(row_dD_px(row_sel)) * p2u;
fprintf('%s :  행 %d 개,  dD 중앙값 %+.3f um  (%s)\n', param.condition, ...
    nnz(row_sel), dD_um, string(unique(row_pol(row_sel))'));

%% 4. A 를 정한다
%  벽이 dD/2 만큼 움직이며 쓸어낸 면적이 2*pi*A 다.
%      2*pi*A = pi*[(r0 + dD/2)^2 - r0^2]
%      A      =    [(r0 + dD/2)^2 - r0^2] / 2
%
%  dD 를 부호째 넣는 것이 핵심이다. 수축이면 dD < 0 이라 (r0 - |dD|/2)^2 이 되고
%  A 가 저절로 음수다. sign(dD) 와 abs(dD) 로 쪼개 쓰면 수축인데 팽창의 기하를
%  적용하게 되어 |A| 가 커진다. see CLAUDE_LOG.md
A_um2 = ((r0_um + dD_um/2)^2 - r0_um^2) / 2;
A_px2 = A_um2 / p2u^2;
flux_true_um2 = 2*pi*A_um2;

fprintf('A          %+.4f um^2\n', A_um2);
fprintf('2 pi A     %+.4f um^2    벽이 쓸어낸 면적. 어느 반경에서도 이만큼 지나간다\n', ...
    flux_true_um2);
fprintf('|u| at 20 um = %.4f um = %.3f px\n', abs(A_um2)/20, abs(A_um2)/20/p2u);

%% 5. 워프를 얹을 프레임 — 조용한 창
%  밑에서 아무것도 안 움직이는 구간이어야 팬텀이 유일한 운동이 된다
quiet_k = find(strcmp({event_det.eventlist.pol}, 'none'), 1);
src_frames = event_det.eventlist(quiet_k).from + (0:param.n_pair-1);
fprintf('조용한 이벤트 %d 번,  프레임 %s\n', quiet_k, mat2str(src_frames));

%% 6. 팬텀 프레임 쌍을 만든다
%  TO 만 보간을 거치므로 FROM 보다 매끄럽다. 그 비대칭은 이득을 낮추는 쪽으로
%  작용한다. see CLAUDE_LOG.md
[grid_x, grid_y] = meshgrid(1:W, 1:H);
off_x = grid_x - centroid_x;
off_y = grid_y - centroid_y;
rr2 = off_x.^2 + off_y.^2;
rr2(rr2 < 1) = 1;
warp_u = A_px2 .* off_x ./ rr2;
warp_v = A_px2 .* off_y ./ rr2;

inter_synth = zeros(H, W, 2*param.n_pair, 'single');
for k = 1:param.n_pair
    frame_from = double(S(:, :, src_frames(k)));
    frame_to = interp2(grid_x, grid_y, frame_from, grid_x - warp_u, grid_y - warp_v, ...
        'cubic', NaN);
    frame_to(isnan(frame_to)) = frame_from(isnan(frame_to));
    inter_synth(:, :, 2*k-1) = single(frame_from);
    inter_synth(:, :, 2*k)   = single(frame_to);
end
fprintf('팬텀 스택 %s.  홀수 = FROM (생), 짝수 = TO (워프)\n', mat2str(size(inter_synth)));
fprintf('워프 크기  최대 %.3f px,  최소 %.3f px\n', ...
    max(hypot(warp_u(:), warp_v(:))), min(hypot(warp_u(:), warp_v(:))));

%% 7. 확인 — 만든 장이 정말 발산 0 인가
%  해석적으로는 정확히 0 이다. 격자에서 재면 유한차분의 절단오차만 남아야 한다
div_num = divergence(grid_x, grid_y, warp_u, warp_v);
outside = ~imdilate(coremask, strel('disk', 8));
fprintf('격자 위 |divergence| 중앙값 (마스크 밖) : %.3e\n', median(abs(div_num(outside))));

figure('Color', 'w', 'Name', 'phantom field');
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile;
step_q = 14;
sub_r = 1:step_q:H;
sub_c = 1:step_q:W;
hold on
visboundaries(coremask, 'Color', [0.85 0.33 0.15], 'LineWidth', 1.4, ...
    'EnhanceVisibility', false);
quiver(grid_x(sub_r, sub_c), grid_y(sub_r, sub_c), ...
    warp_u(sub_r, sub_c)*25, warp_v(sub_r, sub_c)*25, 0, 'Color', [0.1 0.4 0.85]);
axis image
set(gca, 'YDir', 'reverse');
title(sprintf('u = A/r,  A = %+.3f um^2   (x25)', A_um2));
hold off
nexttile;
r_line = linspace(r_wall_um, 45, 300);
plot(r_line, 2*pi*abs(A_um2)*ones(size(r_line)), 'LineWidth', 2.2);
ylim([0 2*pi*abs(A_um2)*1.3]);
grid on
xlabel('반경 r  (um)');
ylabel('\Phi(r) = 2\pi r |u_r|   (um^2)');
title('참값 : 어느 원을 지나든 같은 면적');

%% 8. PIV — 파이프라인이 무엇을 돌려주나
piv_corr = piv_corr_ensemble(double(inter_synth), 'window_sizes', param.window_sizes, ...
    'repeat', 1, 'do_pad', 1, 'exclmask', exclmask, 'use_gpu', param.use_gpu);
[u_valid, v_valid] = piv_validate(piv_corr.utable, piv_corr.vtable, piv_corr.corr);
xyuv_all = cat(3, piv_corr.xtable, piv_corr.ytable, u_valid, v_valid);
keep_vec = (piv_corr.typevector == 1) & ~isnan(u_valid) & ~isnan(v_valid);
xyuv_synth = piv_blank(xyuv_all, keep_vec);
fprintf('벡터 %d 개 중 %d 개 살아남음\n', numel(keep_vec), nnz(keep_vec));

%% 9. 이득 k(r) — 되찾은 크기 / 참 크기
%  균일한 이득은 발산을 못 만든다. 반경에 딸린 이득이 만든다. 그래서 이 줄이
%  평평한지가 셀 10 을 결정한다
mag_meas_um = hypot(xyuv_synth(:,:,3), xyuv_synth(:,:,4)) * p2u;
r_vec_um = hypot(xyuv_synth(:,:,1) - centroid_x, xyuv_synth(:,:,2) - centroid_y) * p2u;
ring_edges = param.ring_edges_um;
gain_r_um = nan(1, numel(ring_edges)-1);
gain_k    = nan(1, numel(ring_edges)-1);
for e = 1:numel(ring_edges)-1
    in_ring = r_vec_um >= ring_edges(e) & r_vec_um < ring_edges(e+1) & isfinite(mag_meas_um);
    gain_r_um(e) = mean(ring_edges(e:e+1));
    if nnz(in_ring) < 20
        continue
    end
    mag_ring = median(mag_meas_um(in_ring));
    mag_true = median(abs(A_um2) ./ r_vec_um(in_ring));
    gain_k(e) = mag_ring / mag_true;
end
fprintf('r  um  ');
fprintf('%8.0f', gain_r_um);
fprintf('\nk      ');
fprintf('%8.3f', gain_k);
fprintf('\n평평하면 volume_out 이 0 이다. 오르내리는 폭이 곧 위양성의 크기\n');

%% 10. volume_out — 참값은 웨지마다 0
%  accumulate 는 bin 을 반경으로 누적하므로, 한 웨지의 값은 그 영역 경계를 지나는
%  flux 와 같다 (선형 보간의 발산정리). 내부 삼각형의 값은 상쇄되어 남지 않는다
vp = vfield_polar(xyuv_synth, coremask, p2u, n_wedge = param.n_wedge, ...
    bin_edges_um = param.bin_edges_um, exclmask = exclmask, gated = true);
vp.param.min_tri_wedge = param.min_tri_wedge;
vp.param.verbose = false;
vp.applydelaunay();
vp.trifilter();
vp.placetri();
vp.measure();
cells = vp.accumulate(vp.gate_wedge());

wedge_vol = cells.volume_out(param.bin_read, :);
fprintf('\n웨지별 volume_out at bin %d   (참값은 웨지마다 0)\n', param.bin_read);
fprintf('   웨지 ');
fprintf('%9d', 1:param.n_wedge);
fprintf('\n   값   ');
fprintf('%+9.4f', wedge_vol);
fprintf('\n   평균 %+.4f   sd %.4f\n', mean(wedge_vol, 'omitnan'), std(wedge_vol, 'omitnan'));
fprintf('\n0 에서 벗어난 만큼이 전부 파이프라인이 만든 것이다\n');
