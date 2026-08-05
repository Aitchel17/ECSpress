function [ax, h] = piv_review_events(events, det, dtrace, t_axis, opt)
%PIV_REVIEW_EVENTS  Review figure for piv_merge_bouts -> piv_pick_excursions output.
%   The gate before PIV: a mis-picked endpoint shows up here.
%   Raw trace plus dsg, the trace the picking saw and the markers sit on. Bands are
%   the raw 5-category scoring in SECONDS; dashed lines are the merged bout edges
%   detection cut on, in FRAMES.
%
% IN   events, det      outputs of piv_pick_excursions
%      dtrace, t_axis   1xT each
%      sleep_score   raw struct for the bands ([] = none)
%      state_idx     frame-index bout tables for the edge lines ([] = none)
%      states        which state_idx fields form the merge partition
%      band_alpha    0.40, as analysis_main
%      ax            target axes ([] = new figure)
% OUT  ax, h (raw, dsg, bands, edges)

arguments
    events  struct
    det     struct
    dtrace  (1,:) double
    t_axis  (1,:) double
    opt.sleep_score = []
    opt.state_idx   = []
    opt.states      cell = {'roughnrem', 'rem', 'roughawake'}
    opt.band_alpha  (1,1) double = 0.40
    opt.ax          = []
end

% 0. Setup
nT = numel(dtrace);
if numel(det.dsg) < nT
    error('piv_review_events:traceShort', ...
        'det.dsg is shorter than dtrace (%d).', nT);
end
if isempty(opt.ax); figure('Name', 'event review'); ax = gca; else; ax = opt.ax; end

% 1. Draw the raw and picking traces
h.raw = plot(ax, t_axis(1:nT), dtrace, 'Color', [.72 .72 .72], 'LineWidth', 0.5);
hold(ax, 'on'); grid(ax, 'on');
h.dsg = plot(ax, t_axis(1:nT), det.dsg(1:nT), 'Color', [.85 0 .85], 'LineWidth', 1.6);
xlim(ax, t_axis([1 nT]));       % sleep_score is untrimmed; keep it off the axis limits
yl = ylim(ax);

% 2. Sleep state bands
% 2.1 sleep_score is in seconds, so no t_axis argument. Colour comes from the field name
h.bands = gobjects(0);
if ~isempty(opt.sleep_score)
    h.bands = plot_sleep_patches(ax, opt.sleep_score, 'YLim', yl, 'FaceAlpha', opt.band_alpha);
    uistack(h.bands, 'bottom');
end

% 3. Merged bout edges
% 3.1 Collect the edges of the merge partition states
h.edges = gobjects(0);
if ~isempty(opt.state_idx)
    medge = [];
    for s = 1:numel(opt.states)
        if isfield(opt.state_idx, opt.states{s})
            medge = [medge; opt.state_idx.(opt.states{s})]; %#ok<AGROW>
        end
    end
    % 3.2 One bout's end is the next one's start, so unique() drops the duplicates
    medge = unique(medge(:))';
    medge = medge(medge >= 1 & medge <= nT);
    % 3.3 Dashed line per edge, frame indices -> seconds
    for x = medge
        h.edges(end+1) = xline(ax, t_axis(x), '--', 'Color', [.35 .35 .35], ...
            'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end

% 4. Event segments: red up = dilation, blue down = constriction
mtr   = det.dsg;
isdil = strcmp({events.pol}, 'dilation');
for e = 1:numel(events)
    if isdil(e); sty = 'r-^'; fc = 'r'; else; sty = 'b-v'; fc = 'b'; end
    plot(ax, t_axis([events(e).from events(e).to]), mtr([events(e).from events(e).to]), ...
        sty, 'LineWidth', 1.6, 'MarkerFaceColor', fc, 'MarkerIndices', 2, ...
        'HandleVisibility', 'off');
end

% 5. Labels
ylim(ax, yl); xlabel(ax, 't (s)'); ylabel(ax, 'diameter (\mum)');
legend(ax, [h.raw h.dsg], ...
    {'raw', sprintf('sgolay %gs (picking)', det.P.sg_win_s)}, 'Location', 'best');
title(ax, sprintf(['%d events: %d dilation (red) / %d constriction (blue)   |   ' ...
    'dashed = merged bout edges'], numel(events), nnz(isdil), nnz(~isdil)));
end
