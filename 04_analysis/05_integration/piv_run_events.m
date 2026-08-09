function [piv_run, info] = piv_run_events(events, S, opt)
%PIV_RUN_EVENTS  Correlation-ensemble PIV for each detected diameter excursion.
%   One displacement field per event, between the two quasi-static states it links.
%   Framesets from +/- halfwin are interleaved as (from_1,to_1,from_2,to_2,...), so
%   every pair spans the whole transition and uv is the total from->to displacement
%   already -- do NOT multiply by nfr.
%
% IN   events   struct array from fwhm_pick_excursions (from/to, plus
%               state/pol/rise_s/diameter_change for the record)
%      S        H x W x T preprocessed stack (piv_preprocess)
%      halfwin  frames either side of from/to forming each frameset (default 2)
%      sel      indices into events ([] = all)
%      piv_opt  cell of name-values passed to piv_corr_ensemble; exclmask belongs here
%      save_planes  keep the per-pair correlation planes (large; default false)
%      verbose  progress line per event (default true)
% OUT  piv_run  1xN struct: id, state, pol, from, to, rise_s, diameter_change,
%               nfr (pairs averaged, not a scale factor), uv (HxWx2, px), corr,
%               planes
%      info     halfwin, piv_opt, sel, nfr

arguments
    events struct
    S      {mustBeNumeric, mustBeNonempty}
    opt.halfwin     (1,1) double  = 2
    opt.sel         double  = []
    opt.piv_opt     cell    = {'window_sizes', [40 20; 20 10; 12 6], 'repeat', 1, 'do_pad', 1}
    opt.save_planes (1,1) logical = false
    opt.verbose     (1,1) logical = true
end

% 0. Setup
% 0.1 Event selection
nT  = size(S, 3);
sel = opt.sel;
if isempty(sel); sel = 1:numel(events); end
sel = reshape(sel, 1, []);
if any(sel < 1 | sel > numel(events))
    error('piv_run_events:badSel', 'sel indexes outside events (1..%d).', numel(events));
end
% 0.2 PIV options
popt = opt.piv_opt;
if opt.save_planes; popt = [popt, {'save_corrmaps', true}]; end
% 0.3 Output template
piv_run = struct('id',{},'state',{},'pol',{},'from',{},'to',{},'rise_s',{}, ...
                 'diameter_change',{},'nfr',{},'uv',{},'corr',{},'planes',{});

% 1. One ensemble PIV run per selected event
for n = 1:numel(sel)
    e = sel(n);
    E = events(e);
    % 1.1 The event must run forward in time (pair order sets the sign) and the two
    %     framesets must not overlap, or pivensemble correlates a frame with itself
    if E.to - E.from <= 2 * opt.halfwin
        error('piv_run_events:framesetOverlap', ...
            ['event %d spans %d frames (from %d to %d); halfwin=%d needs more ' ...
             'than %d.'], e, E.to - E.from, E.from, E.to, opt.halfwin, 2 * opt.halfwin);
    end
    % 1.2 Framesets around each endpoint
    w_from = max(1, E.from - opt.halfwin):min(nT, E.from + opt.halfwin);
    w_to   = max(1, E.to   - opt.halfwin):min(nT, E.to   + opt.halfwin);
    if opt.verbose
        fprintf('[%d/%d] event %d  %-14s %-12s %5.1f s ', n, numel(sel), e, ...
            E.state, E.pol, E.rise_s);
    end
    % 1.3 Interleave and run
    [uv, cc, cp] = ens_interleave(S, w_from, w_to, popt);
    piv_run(n) = struct('id', e, 'state', E.state, 'pol', E.pol, ...
        'from', E.from, 'to', E.to, 'rise_s', E.rise_s, ...
        'diameter_change', E.diameter_change, ...
        'nfr', min(numel(w_from), numel(w_to)), 'uv', uv, 'corr', cc, 'planes', cp);
    if opt.verbose
        fprintf('| %+.2f um, %d pairs\n', E.diameter_change, piv_run(n).nfr);
    end
end

% 2. Settings record
info = struct('halfwin', opt.halfwin, 'piv_opt', {popt}, 'sel', sel, ...
              'nfr', [piv_run.nfr]);
end

% ---------------------------------------------------------------------------
function [uv, cc, cp] = ens_interleave(S, w_from, w_to, popt)
%ENS_INTERLEAVE  Interleave the earlier and later framesets -> correlation ensemble.
%   w_from must be the earlier frameset: the pair order sets the sign of uv. Pair p
%   is (from_p, to_p), so the ensemble averages k measurements of one displacement
%   rather than accumulating k steps.
    % 0. Equal-length index lists for the two framesets
    k = min(numel(w_from), numel(w_to));
    a = round(linspace(w_from(1), w_from(end), k));
    b = round(linspace(w_to(1),   w_to(end),   k));
    % 1. Interleave into (from_1, to_1, from_2, to_2, ...)
    inter = zeros(size(S, 1), size(S, 2), 2 * k);
    inter(:, :, 1:2:end) = S(:, :, a);
    inter(:, :, 2:2:end) = S(:, :, b);
    % 2. Ensemble PIV. It correlates only, so the gating is explicit here
    cp = piv_corr_ensemble(inter, popt{:});
    [cp.utable, cp.vtable] = piv_validate(cp.utable, cp.vtable, cp.corr);
    [uv, cc] = vfield_stamp(cp);
    % 3. Relabel pairs with original stack indices
    cp.pair_frames = [a(:), b(:)];
end
