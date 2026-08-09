function [ax, h] = event_review(events, dtrace, smooth_diameter, sg_win_s, t_axis, opt)
%EVENT_REVIEW  Review figure for event_merge_bouts -> event_pick_excursions output.
%   The gate before PIV: a mis-picked endpoint shows up here.
%   Raw trace plus dsg, the trace the picking saw and the markers sit on. Bands are
%   the raw 5-category scoring in SECONDS; dashed lines are the merged bout edges
%   detection cut on, in FRAMES.
%
% IN   events           1 x N struct  event_pick_excursions output, controls included
%      dtrace           1 x T float   raw caliber trace
%      smooth_diameter  1 x T float   the trace the picking saw; markers sit on it
%      sg_win_s         float s       Savitzky-Golay window, for the legend only
%      t_axis           1 x T float   frame times (s)
%      pixel2um      float um/px on both traces, default 1 = they are
%                    already in um. This is the figure layer, so the conversion
%                    happens here -- same place analysis_pax_makefig does it
%      sleep_score   raw struct for the bands ([] = none)
%      state_idx     frame-index bout tables for the edge lines ([] = none)
%      states        which state_idx fields form the merge partition
%      focus         cell of event labels to draw in POLARITY colour; every other
%                    transition goes thin grey. {} = colour them all
%      show_control  false hides the pol='none' stable windows entirely
%      band_alpha    0.40, as analysis_main
%      ax            target axes ([] = new figure)
% OUT  ax, h (raw, dsg, bands, edges)
%
%   why  focus greys instead of dropping : a picked endpoint is judged against the
%        bouts AROUND it, so the neighbouring transitions have to stay on the axes
%        even when they are not the ones being read

arguments
    events  struct
    dtrace           (1,:) double
    smooth_diameter  (1,:) double
    sg_win_s         (1,1) double
    t_axis  (1,:) double
    opt.pixel2um    (1,1) double = 1
    opt.sleep_score = []
    opt.state_idx   = []
    opt.states      cell = {'roughnrem', 'rem', 'roughawake'}
    opt.focus       cell = {}
    opt.show_control (1,1) logical = true
    opt.band_alpha  (1,1) double = 0.40
    opt.ax          = []
end
dtrace          = dtrace          * opt.pixel2um;
smooth_diameter = smooth_diameter * opt.pixel2um;

% 0. Setup
nT = numel(dtrace);
if numel(smooth_diameter) < nT
    error('event_review:traceShort', ...
        'smooth_diameter is shorter than dtrace (%d).', nT);
end
if isempty(opt.ax); figure('Name', 'event review'); ax = gca; else; ax = opt.ax; end

% 1. Draw the raw and picking traces
h.raw = plot(ax, t_axis(1:nT), dtrace, 'Color', [.72 .72 .72], 'LineWidth', 0.5);
hold(ax, 'on'); grid(ax, 'off');
h.dsg = plot(ax, t_axis(1:nT), smooth_diameter(1:nT), 'Color', [.85 0 .85], 'LineWidth', 1.6);
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

% 4. Event segments: red up = dilation, blue down = constriction, grey flat = none
%    err  a two-way isdil test drew every pol that was not 'dilation' as a blue
%         constriction, so the stable windows quietdetection adds came out
%         labelled as the thing they are the control FOR
mtr   = smooth_diameter;                       % already scaled at the top of the file
isdil = strcmp({events.pol}, 'dilation');
iscon = strcmp({events.pol}, 'constriction');
isctl = ~isdil & ~iscon;
if isempty(opt.focus)
    infocus = true(1, numel(events));
else
    infocus = ismember({events.state}, opt.focus);
end
n_drawn = 0;  n_grey = 0;
for e = 1:numel(events)
    if isctl(e) && ~opt.show_control
        continue
    end
    % err  the out-of-focus grey started at [.62 .62 .62] and vanished : the awake
    %      sleep band is grey at the same lightness, so 28 segments drew invisibly
    %      on top of it. Near-black reads over every band colour and still sits
    %      clearly behind the red/blue without competing
    if isctl(e);           sty = '-s';  fc = [.45 .45 .45];
    elseif ~infocus(e);    sty = '-';   fc = [.15 .15 .15];
    elseif isdil(e);       sty = 'r-^'; fc = 'r';
    else;                  sty = 'b-v'; fc = 'b';
    end
    hseg = plot(ax, t_axis([events(e).from events(e).to]), ...
        mtr([events(e).from events(e).to]), sty, 'LineWidth', 1.6, ...
        'MarkerFaceColor', fc, 'MarkerIndices', 2, 'HandleVisibility', 'off');
    if isctl(e)
        set(hseg, 'Color', fc, 'LineWidth', 2.4);
        n_drawn = n_drawn + 1;
    elseif ~infocus(e)
        set(hseg, 'Color', fc, 'LineWidth', 1.2);
        n_grey = n_grey + 1;
    else
        n_drawn = n_drawn + 1;
    end
end

% 5. Labels
ylim(ax, yl); xlabel(ax, 't (s)'); ylabel(ax, 'diameter (\mum)');
legend(ax, [h.raw h.dsg], ...
    {'raw', sprintf('sgolay %gs (picking)', sg_win_s)}, 'Location', 'best');
% the title states what was HIDDEN as well as what was drawn -- a review figure
% that quietly dropped rows would be the one thing this plot exists to prevent
parts = sprintf('%d dilation (red) / %d constriction (blue)', ...
    nnz(isdil & infocus), nnz(iscon & infocus));
if n_grey > 0;  parts = sprintf('%s / %d other (grey)', parts, n_grey); end
if opt.show_control
    parts = sprintf('%s / %d stable (grey square)', parts, nnz(isctl));
elseif any(isctl)
    parts = sprintf('%s   |   %d stable HIDDEN', parts, nnz(isctl));
end
title(ax, sprintf('%d of %d rows: %s   |   dashed = merged bout edges', ...
    n_drawn + n_grey, numel(events), parts));
end
