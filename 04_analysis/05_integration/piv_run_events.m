function [piv_run, info] = piv_run_events(events, stack, fps, pixel2um, opt)
%PIV_RUN_EVENTS  One analysis_pivensemble per detected event : cut - correlate - gate - row
%   The frames are cut by piv_getframescope, the rest is the class. uv is the total
%   from->to displacement already -- do NOT multiply by nfr.
%
% IN   events    struct array from analysis_event (from / to, plus state / pol / bout /
%                rise_s / diameter_change for the record)
%      stack     H x W x T numeric  the recording before PIV preprocessing
%      fps       1 x 1 double
%      pixel2um  1 x 1 double
%      halfwin   1 x 1 int       frames either side of each endpoint (default 2)
%      sel       1 x n int       indices into events ([] = all)
%      windows   P x 2 int       [window step] per pass, into param.pivlab_windows
%      exclmask  H x W bool      true = pixel excluded, into param.pivlab_mask ([] = none)
%      verbose   1 x 1 bool      one line per event (default true)
% OUT  piv_run   1 x N struct : id, state, pol, bout, from, to, rise_s, diameter_change,
%                nfr (pairs averaged, not a scale factor), pair_frames (nfr x 2, recording
%                frames), xyuv (ny x nx x 4 px, gated), xyuv_ungated, common_mode (1 x 2 px),
%                gate_name, gates, planes (grid, vectors, corr, typevector, peak_floor; no maps)
%      info      halfwin, windows, sel, nfr, gate_name

arguments
    events struct
    stack  {mustBeNumeric, mustBeNonempty}
    fps      (1,1) double {mustBePositive}
    pixel2um (1,1) double {mustBePositive}
    opt.halfwin  (1,1) double = 2
    opt.sel      double = []
    opt.windows  (:,2) double = [40 20; 20 10; 12 6]
    opt.exclmask logical = []
    opt.verbose  (1,1) logical = true
end

% 0. event selection
sel = opt.sel;
if isempty(sel)
    sel = 1:numel(events);
end
sel = reshape(sel, 1, []);
if any(sel < 1 | sel > numel(events))
    error('piv_run_events:badSel', 'sel indexes outside events (1..%d)', numel(events));
end
% bout is carried : two transitions out of one bout are not independent, and that
% is unrecoverable once the row is written
piv_run = struct('id', {}, 'state', {}, 'pol', {}, 'bout', {}, 'from', {}, 'to', {}, 'rise_s', {}, ...
                 'diameter_change', {}, 'nfr', {}, 'pair_frames', {}, 'xyuv', {}, 'xyuv_ungated', {}, ...
                 'common_mode', {}, 'gate_name', {}, 'gates', {}, 'planes', {});

% 1. one object per selected event
for n = 1:numel(sel)
    e = sel(n);
    E = events(e);
    if opt.verbose
        fprintf('[%d/%d] event %d  %-14s %-12s %5.1f s ', n, numel(sel), e, E.state, E.pol, E.rise_s);
    end
    [stack_span, span] = piv_getframescope(stack, E.from, E.to, opt.halfwin);
    ensemble = analysis_pivensemble(stack_span, opt.halfwin, fps, pixel2um);
    ensemble.param.pivlab_windows = opt.windows;
    ensemble.param.pivlab_mask    = opt.exclmask;
    ensemble.param.verbose        = false;
    ensemble.correlate("endpoint");
    ensemble.gate("endpoint");
    endpoint = ensemble.endpoint;
    piv_run(n) = struct('id', e, 'state', E.state, 'pol', E.pol, ...
        'bout', E.bout, 'from', E.from, 'to', E.to, 'rise_s', E.rise_s, ...
        'diameter_change', E.diameter_change, ...
        'nfr', size(endpoint.pair_frames, 1), ...
        'pair_frames', span(endpoint.pair_frames), ...
        'xyuv', endpoint.xyuv, 'xyuv_ungated', endpoint.xyuv_ungated, ...
        'common_mode', endpoint.common_mode, 'gate_name', endpoint.gate_name, ...
        'gates', endpoint.gates, 'planes', endpoint.planes);
    if opt.verbose
        fprintf('| %+.2f um, %d pairs, %d vectors\n', E.diameter_change, piv_run(n).nfr, ...
            nnz(~isnan(endpoint.xyuv(:,:,3))));
    end
end

% 2. settings record
info = struct('halfwin', opt.halfwin, 'windows', opt.windows, 'sel', sel, ...
              'nfr', [piv_run.nfr], 'gate_name', string({piv_run.gate_name}));
end
