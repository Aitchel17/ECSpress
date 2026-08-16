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
% OUT  piv_run  1xN struct: id, state, pol, bout, from, to, rise_s, diameter_change,
%               nfr (pairs averaged, not a scale factor),
%               xyuv (ny x nx x 4 px, [x y u v] on the interrogation grid),
%               xyuv_ungated (the same vectors before piv_validate),
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
% bout says which sleep bout the event sits in. Carried, not dropped: two
% transitions of the same kind out of the SAME bout are less independent than
% two out of different ones, and that is unrecoverable once the row is written
piv_run = struct('id',{},'state',{},'pol',{},'bout',{},'from',{},'to',{},'rise_s',{}, ...
                 'diameter_change',{},'nfr',{},'xyuv',{},'xyuv_ungated',{},'planes',{});

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
    if opt.verbose
        fprintf('[%d/%d] event %d  %-14s %-12s %5.1f s ', n, numel(sel), e, ...
            E.state, E.pol, E.rise_s);
    end
    % 1.2 Interleave and run
    [xyuv, xyuv_ungated, cp] = ens_interleave(S, E.from, E.to, opt.halfwin, popt);
    piv_run(n) = struct('id', e, 'state', E.state, 'pol', E.pol, ...
        'bout', E.bout, 'from', E.from, 'to', E.to, 'rise_s', E.rise_s, ...
        'diameter_change', E.diameter_change, ...
        'nfr', size(cp.pair_frames, 1), 'xyuv', xyuv, ...
        'xyuv_ungated', xyuv_ungated, 'planes', cp);
    if opt.verbose
        fprintf('| %+.2f um, %d pairs\n', E.diameter_change, piv_run(n).nfr);
    end
end

% 2. Settings record
info = struct('halfwin', opt.halfwin, 'piv_opt', {popt}, 'sel', sel, ...
              'nfr', [piv_run.nfr]);
end

% ---------------------------------------------------------------------------
function [xyuv, xyuv_ungated, cp] = ens_interleave(S, from, to, halfwin, popt)
%ENS_INTERLEAVE  Correlate the pairs piv_interleave makes, then read the verdict.
%   The interleave itself is piv_interleave, called and not copied: a second copy
%   of that arithmetic would decide which frames a result came from, and drift in
%   it would make "the frames PIV used" a false statement with nothing to catch it.
    % 0. The pairs
    [inter, pair_frames] = piv_interleave(S, from, to, halfwin);
    % 1. Ensemble PIV. It correlates only, so the gating is explicit here
    cp = piv_corr_ensemble(inter, popt{:});
    % piv_validate answers by NaN-ing what it rejects, so the raw field is kept
    % first and the verdict read back off the difference. Both maps then come off
    % the SAME vectors and differ only in their mask
    u_raw = cp.utable;
    v_raw = cp.vtable;
    [cp.utable, cp.vtable] = piv_validate(u_raw, v_raw, cp.corr);
    xyuv_all   = cat(3, cp.xtable, cp.ytable, u_raw, v_raw);
    keep_raw   = (cp.typevector == 1) & ~isnan(u_raw) & ~isnan(v_raw);
    keep_gated = keep_raw & ~isnan(cp.utable) & ~isnan(cp.vtable);
    % the row carries the grid. The dense map is a view and piv_stamp makes it
    % where one is wanted, which is the display path only
    xyuv         = piv_blank(xyuv_all, keep_gated);
    xyuv_ungated = piv_blank(xyuv_all, keep_raw);
    % 3. Relabel pairs with original stack indices
    cp.pair_frames = pair_frames;
end

