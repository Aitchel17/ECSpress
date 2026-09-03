%PIV_VALIDATION_TESTBED  Why the correlation plane cannot decide which vector to keep.
%
% The question this answers: the correlation plane is right there, so why gate on
% the IMAGE (Shi-Tomasi lambda2) instead of on the correlation? Every panel is
% see FINDINGS.md
%
% ANSWER, in three parts
%   peak HEIGHT  corr2 falls with displacement, so it spends its budget on the
%                strongest signal. At the same budget as lambda2 it rejects big
%                labelled-GOOD vectors where lambda2 takes none. see FINDINGS.md
%   peak SHAPE   the width ratio is at chance. A textureless window still gives a
%                sharp peak -- confident, and in the wrong place. see FINDINGS.md
%   plane FLOOR  the one plane measure that carries something: the rms AWAY from
%                the peak. Kept as gate_corr_floor, second to lambda2, not instead
%
% caution  "the gate rejects long vectors" is NOT by itself the charge : both gates
%          reject more at large |d| (panel B), because the labelled-BAD windows are
%          the longer ones. The charge is WHICH long ones, panel C. see FINDINGS.md
%
% needs   piv_ensemble, an analysis_pivensemble correlated and gated on ev2
% writes  piv_validation_corrgates.svg next to this file
% caution R2023b exportgraphics has no SVG writer; print -dsvg -vector does
% caution the labels are window groups drawn by eye on ev2's quiver plot. A second
%         event of opposite polarity (ev22) reproduced every conclusion but is not
%         re-run here. see FINDINGS.md

%% 0. Labels. 20 px circles around each centre, drawn by eye on the quiver plot
label_live = {[158 98], [115 85], [80 194]};
label_die  = {[338 134], [254 146], [212 20], [296 200], [68 218], [62 218], ...
        [296 104], [332 98], [362 206], [236 32], [308 26], [50 188]};

result = piv_ensemble.endpoint;
planes = result.planes;
maps   = piv_ensemble.get_corrplane("endpoint");   % the object keeps no planes
mask_ok = ~isnan(planes.utable) & ~isnan(planes.vtable) & planes.typevector == 1;
len_grid = hypot(planes.utable, planes.vtable);
% note  lambda2(:,:,1) = leading frames, (:,:,2) = trailing; the worse end decides
lambda2 = min(result.gates(1).info.lambda2, [], 3);

circle = @(p) mask_ok & hypot(planes.xtable-p(1), planes.ytable-p(2)) <= 20;
cell_live = cellfun(circle, label_live, 'uni', 0);
cell_die  = cellfun(circle, label_die,  'uni', 0);
mask_live = any(cat(3, cell_live{:}), 3);
mask_die  = any(cat(3, cell_die{:}),  3);

%% 1. Every candidate, read off the final-pass plane one window at a time
[ny, nx] = size(planes.utable);
n_cell   = size(maps, 1);              % plane side, px. 1 cell = 1 px
[YY, XX] = ndgrid(1:n_cell, 1:n_cell);
plane_measure = struct( ...
    'peak_mean',  NaN(ny,nx), ...   % ny x nx float  peak over the plane mean
    'peak_mass',  NaN(ny,nx), ...   % ny x nx float  fraction of energy within 2 px
    'corr_floor', NaN(ny,nx), ...   % ny x nx float  peak over rms away from it
    'w_ratio',    NaN(ny,nx));      % ny x nx float  width ALONG d over ACROSS d
for r = 1:ny
    for c = 1:nx
        if ~mask_ok(r,c); continue; end
        plane_E = double(sum(maps(:,:,r,c,:), 5));
        plane_E = plane_E - min(plane_E(:));
        [peak, k] = max(plane_E(:));
        if peak <= 0; continue; end
        [peak_row, peak_col] = ind2sub(size(plane_E), k);
        mask_near = hypot(XX-peak_col, YY-peak_row) <= 2;
        plane_measure.peak_mean(r,c)  = peak / mean(plane_E(:));
        plane_measure.peak_mass(r,c)  = sum(plane_E(mask_near)) / sum(plane_E(:));
        plane_measure.corr_floor(r,c) = peak / max(rms(plane_E(~mask_near)), eps);
        % second moment of the half-max core, resolved onto the vector's own
        % direction. + 1/12 is the variance of a 1 px cell, so a single-cell peak
        % gives ratio 1 instead of 0/0
        core = max(plane_E - 0.5*peak, 0);
        if sum(core(:)) > 0
            w  = core / sum(core(:));
            dx = XX - peak_col;   dy = YY - peak_row;
            M  = [sum(w(:).*dx(:).^2)        sum(w(:).*dx(:).*dy(:)); ...
                  sum(w(:).*dx(:).*dy(:))    sum(w(:).*dy(:).^2)];
            th = atan2(planes.vtable(r,c), planes.utable(r,c));
            e_along  = [cos(th); sin(th)];
            e_across = [-sin(th); cos(th)];
            plane_measure.w_ratio(r,c) = sqrt(e_along'*M*e_along + 1/12) / ...
                                         sqrt(e_across'*M*e_across + 1/12);
        end
    end
end

%% 2. Score every candidate against the labels
% note  plain text, not identifiers : these are axis labels and TeX would read an
%       underscore as a subscript
name_measure = ["corr2", "peak/mean", "peak mass", "width ratio", "corr floor", "\lambda_2"];
val_measure  = {planes.corr, plane_measure.peak_mean, plane_measure.peak_mass, ...
                plane_measure.w_ratio, plane_measure.corr_floor, lambda2};
auc = zeros(1, numel(val_measure));
for m = 1:numel(val_measure)
    a = val_measure{m}(mask_live);   a = a(~isnan(a));
    b = val_measure{m}(mask_die);    b = b(~isnan(b));
    % Mann-Whitney AUC, then folded: a measure that separates the wrong way round
    % still separates, and the sign is a convention
    q = mean(reshape(a,[],1) >  reshape(b,1,[]), 'all') + ...
    0.5*mean(reshape(a,[],1) == reshape(b,1,[]), 'all');
    auc(m) = max(q, 1-q);
end

% 2.1 The three verdicts, as the class actually forms them. corr2 is the RETIRED
%     gate_rheight at PIVlab's own 0.5, so the comparison is against a real gate
mask_rej = { mask_ok & planes.corr < 0.5, ...
             mask_ok & result.gates(1).mask, ...
             mask_ok & result.gates(2).mask };
name_gate = ["corr2 < 0.5  (retired)", ...
             sprintf('\\lambda_2 < %.2g x median', piv_ensemble.param.tomasi_medfrac), ...
             sprintf('corr floor < %.2g x median', piv_ensemble.param.corr_floor_medfrac)];

% 2.2 By displacement decile. Both gates rise -- see the caution in the header
edges_d  = prctile(len_grid(mask_ok), 0:10:100);
d_mid    = zeros(1,10);   corr_mid = zeros(1,10);
pct_dec  = zeros(3,10);
for k = 1:10
    in = mask_ok & len_grid >= edges_d(k) & len_grid <= edges_d(k+1);
    d_mid(k)    = median(len_grid(in));
    corr_mid(k) = median(planes.corr(in));
    for g = 1:3;  pct_dec(g,k) = 100*mean(mask_rej{g}(in));  end
end
r_corr_d = corr(len_grid(mask_ok), planes.corr(mask_ok));

% 2.3 By label AND length -- the panel that decides it. Split at the field median
%     so the four groups are disjoint and cover every labelled window
is_big     = len_grid > median(len_grid(mask_ok), 'omitnan');
mask_group = {mask_live & ~is_big, mask_live & is_big, ...
              mask_die  & ~is_big, mask_die  & is_big};
name_group = ["GOOD short", "GOOD long", "BAD short", "BAD long"];
pct_group  = zeros(4, 3);
pct_die    = zeros(1, 3);            % pooled over BAD, i.e. each gate's budget
for g = 1:3
    pct_die(g) = 100*mean(mask_rej{g}(mask_die));
    for q = 1:4;  pct_group(q,g) = 100*mean(mask_rej{g}(mask_group{q}));  end
end

%% 3. The figure
fig = figure('Color', 'w', 'Position', [30 30 1560 1080]);
tl  = tiledlayout(fig, 4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
CL  = color_lee.clist;
co_gate    = [CL.red; CL.blue; CL.teal];
cmap_plane = color_lee().gradient.inferno;

% A. corr2 falls as the displacement grows
nexttile;
scatter(len_grid(mask_ok), planes.corr(mask_ok), 4, CL.lightgray, 'filled'); hold on
plot(d_mid, corr_mid, '-o', 'Color', CL.red, 'LineWidth', 2, 'MarkerFaceColor', CL.red);
xlabel('|d| (px)');  ylabel('corr2');  grid on
title(sprintf('A  corr2 vs displacement,  r = %+.3f', r_corr_d), 'FontWeight', 'normal');

% B. every gate rises with |d|, so this view cannot tell them apart
nexttile;
hold on; grid on
for g = 1:3
    plot(1:10, pct_dec(g,:), '-o', 'Color', co_gate(g,:), 'LineWidth', 1.8, ...
        'MarkerFaceColor', co_gate(g,:), 'MarkerSize', 4);
end
xlabel('|d| decile');  ylabel('rejected (%)');
legend(name_gate, 'Location', 'northwest', 'FontSize', 7);
title('B  all three rise with |d| -- not the discriminator', 'FontWeight', 'normal');

% C. the discriminator. Same budget, different victims
nexttile;
h_bar = bar(pct_group, 'EdgeColor', 'k');
for g = 1:3;  h_bar(g).FaceColor = co_gate(g,:);  end
xticklabels(compose('%s (n=%d)', name_group(:), cellfun(@nnz, mask_group(:))));
ylabel('rejected (%)');  grid on
legend(name_gate, 'Location', 'northwest', 'FontSize', 7);
title(sprintf('C  corr2 takes %.0f%% of the long GOOD vectors, \\lambda_2 %.0f%%', ...
    pct_group(2,1), pct_group(2,2)), 'FontWeight', 'normal');

% D/E. the labels themselves: corr2 overlaps, lambda2 does not
for m = [1 6]
    nexttile;
    a = val_measure{m}(mask_live);   b = val_measure{m}(mask_die);
    histogram(a(~isnan(a)), 20, 'FaceColor', CL.blue, 'FaceAlpha', 0.55, ...
        'Normalization', 'probability');  hold on
    histogram(b(~isnan(b)), 20, 'FaceColor', CL.red,  'FaceAlpha', 0.55, ...
        'Normalization', 'probability');
    xlabel(name_measure(m));  ylabel('fraction');  grid on
    legend({'labelled GOOD', 'labelled BAD'}, 'Location', 'best', 'FontSize', 7);
    title(sprintf('%s  %s,  AUC %.2f', char('D' + find([1 6] == m) - 1), ...
        name_measure(m), auc(m)), 'FontWeight', 'normal');
end

% F. every candidate ranked. Coloured = what the class actually gates on
nexttile;
[auc_sorted, ord] = sort(auc);
barh(auc_sorted, 'FaceColor', CL.lightgray, 'EdgeColor', 'k');  hold on
for k = 1:numel(ord)
    switch name_measure(ord(k))
        case "\lambda_2";   co = CL.blue;
        case "corr floor";  co = CL.teal;
        otherwise;          continue
    end
    barh(k, auc_sorted(k), 'FaceColor', co, 'EdgeColor', 'k');
end
xline(0.5, 'k--');
yticks(1:numel(auc));  yticklabels(name_measure(ord));  xlim([0.4 1]);
xlabel('AUC, labelled GOOD vs BAD');  grid on
title('F  coloured = kept as a gate, grey = rejected', 'FontWeight', 'normal');

% G. the conclusion in text, so the figure stands alone
%    caution  the tile is half the figure wide and a quarter tall. Lines past ~14
%             at this size run out under the plane row
nexttile([1 2]);
axis off
text(0, 1, { ...
 '\bfG   Peak HEIGHT spends its budget on the signal.\rm', ...
 sprintf('     corr2 falls with |d| (r = %+.3f) -- a window that moved far overlaps itself less. Both gates reject more at large |d| (B),', r_corr_d), ...
 sprintf('     and so they should: the labelled-BAD windows ARE the longer ones (median %.2f px against %.2f). What separates the two', ...
          median(len_grid(mask_die)), median(len_grid(mask_live))), ...
 sprintf('     is WHICH long ones. At the same budget (%.0f%% and %.0f%% of labelled-BAD) corr2 takes \\bf%.0f%%\\rm of the long labelled-GOOD', ...
          pct_die(1), pct_die(2), pct_group(2,1)), ...
 sprintf('     vectors and \\lambda_2 takes \\bf%.0f%%\\rm  (C, n = %d).', pct_group(2,2), nnz(mask_group{2})), '', ...
 '\bf    Peak SHAPE carries nothing at all.\rm', ...
 sprintf('     width ratio AUC %.2f -- chance. peak mass %.2f, peak/mean %.2f, against \\lambda_2 %.2f on the same %d labelled windows.', ...
          auc(4), auc(3), auc(2), auc(6), nnz(mask_live) + nnz(mask_die)), ...
 '     H and I are eight planes at matched corr2 quantiles inside each label group: two rows that look alike, labelled opposite.', '', ...
 '\bf    The question is answered BEFORE the correlation runs.\rm', ...
 '     A textureless window still correlates sharply, against whatever it does contain. The peak is confident and it is in the wrong', ...
 '     place, and nothing in the plane says so. Only the IMAGE can, and Shi-Tomasi \lambda_2 is that measurement.', '', ...
 sprintf('\\bf    One plane measure survives:\\rm the rms AWAY from the peak (AUC %.2f). Weaker than \\lambda_2 and overlapping it by only', auc(5)), ...
 '     22-29%, so it is kept as gate\_corr\_floor -- a second axis, not a replacement.'}, ...
 'VerticalAlignment', 'top', 'FontSize', 7.5);

% H/I. eight planes at matched corr2 quantiles: GOOD on one row, BAD on the next
%      caution  SELECTED, at the 20/40/60/80th corr2 percentile of each group, so
%               the two rows span the same range of the measure under test.
%               Neither the best nor the worst of either group
quantile_pick = [0.2 0.4 0.6 0.8];
for row = 1:2
    if row == 1; lin = find(mask_live); tag = 'H  labelled GOOD';
    else;        lin = find(mask_die);  tag = 'I  labelled BAD';
    end
    [~, ord_corr] = sort(planes.corr(lin));
    pick = lin(ord_corr(max(1, round(numel(ord_corr) * quantile_pick))));
    for q = 1:4
        nexttile;
        [r, c] = ind2sub([ny nx], pick(q));
        plane_E = double(sum(maps(:,:,r,c,:), 5));
        plane_E = plane_E - min(plane_E(:));
        imagesc((1:n_cell)-(n_cell+1)/2, (1:n_cell)-(n_cell+1)/2, plane_E);
        axis image;  colormap(gca, cmap_plane);
        xlabel('u residual (px)');
        if q == 1; ylabel({tag, 'v residual (px)'}, 'FontWeight', 'bold');
        else;      ylabel('v residual (px)');
        end
        title({sprintf('(%.0f,%.0f)  |d| %.2f px', planes.xtable(r,c), ...
                       planes.ytable(r,c), len_grid(r,c)), ...
               sprintf('corr2 %.2f  pk/mean %.0f  \\lambda_2 %.1e', ...
                       planes.corr(r,c), plane_measure.peak_mean(r,c), lambda2(r,c))}, ...
              'FontWeight', 'normal', 'FontSize', 8);
        set(gca, 'FontSize', 7);
    end
end

title(tl, sprintf(['Correlation-plane gating vs the Shi-Tomasi image gate   |   ' ...
    'ev2, %d windows correlated, %d labelled GOOD, %d labelled BAD'], ...
    nnz(mask_ok), nnz(mask_live), nnz(mask_die)), 'FontWeight', 'bold');

%% 4. Out
out = fullfile(fileparts(mfilename('fullpath')), 'piv_validation_corrgates.svg');
print(fig, out, '-dsvg', '-vector');
fprintf('wrote %s\n', out);
fprintf('AUC   %s\n', strjoin(compose("%s %.2f", name_measure(:), auc(:)), ' | '));
fprintf('%-22s %s\n', 'rejected (%)', sprintf('%-11s', name_group));
for g = 1:3
    fprintf('%-22s %s\n', name_gate(g), sprintf('%-11.0f', pct_group(:,g)));
end
