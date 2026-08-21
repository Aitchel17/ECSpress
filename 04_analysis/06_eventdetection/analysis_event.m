classdef analysis_event < handle
%ANALYSIS_EVENT  Arousal-state transitions and their controls, as one event list.
%   Takes a diameter series on its frame axis and the bout tables that say which
%   arousal state the animal was in, and returns the spans PIV should run on.
%
%   A row is (state, pol), two axes and no third category:
%     'nrem2rem', 'dilation'      an arousal-state transition
%     'nrem',     'none'          a stable stretch inside one nrem bout, its control
%     'whisker',  'dilation'      a stimulus response, when that method exists
%
%   THE OBJECT IS THE TRACE, not the experiment. Only the diameter and its frame
%   axis are constructor arguments; what defines a span comes in per method. A
%   whisker-stim recording has no sleep bouts and would take stimulus times
%   instead, and it still fills the same eventlist.
%
%   ORDER MATTERS. control() needs arousaltransition() to have run: the controls
%   share its smoothed trace, and the transitions are what they have to stay clear
%   of. control() restarts from the transitions every call, so calls do not stack.
%
%   PX THROUGHOUT, like every analysis class here; the *_makefig layer converts.
%   The picking is scale-invariant, so this costs nothing -- see
%   event_pick_excursions.

    properties
        diameter        (1,:) double  % 1 x T float   caliber series, PX
        t_axis          (1,:) double  % 1 x T float   frame times (s); fps derives from this
        smooth_diameter (1,:) double  % 1 x T float   what the picking ran on; set by arousaltransition
        eventlist                     % 1 x N struct  every event and control, chronological
        excursion       table         % every diameter excursion, enumerated ONCE.
                                      % see design/EXCURSION_DESIGN.md
        info            struct        % what each stage measured, see below
        param = struct( ...           % the partition. Both stages read it
            'merge', {{ ...           % N x 4 cell  {earlier, later, label, polarity}
        %    earlier(A)    later(B)      label             polarity
            'nrem',       'rem',        'nrem2rem',       'dilation'     ; ...
            'uarousal',   'nrem',       'uarousal2nrem',  'dilation'     ; ...
            'roughawake', 'nrem',       'awake2nrem',     'dilation'     ; ...
            'rem',        'roughawake', 'rem2awake',      'constriction' ; ...
            'nrem',       'roughawake', 'nrem2awake',     'constriction' ; ...
            'nrem',       'uarousal',   'nrem2uarousal',  'constriction' }}, ...
            'states', {{'roughnrem', 'rem', 'roughawake'}})
                                      % 1 x K cell  state_bout fields the controls search
    end
    %   why param and not two arguments : control's own header requires the SAME
    %   partition arousaltransition was given -- it finds boundaries BETWEEN states
    %   and control searches INSIDE them, so a different partition makes the
    %   avoidance meaningless. Passed separately at two call sites, nothing shows
    %   whether they agree. Held here, one object answers it
    %   caution  THIS IS THE ONE PLACE THE CLASS NAMES SLEEP STATES. Everything
    %            else about scoring stays outside -- the bout tables come in as an
    %            argument, and nothing here interprets them. The default is this
    %            project's standard six transitions; a design with different states
    %            sets param.merge and param.states and the class does not care
    %   note  merge and states do not have to name the same fields, and by default
    %         they do not: merge names nrem / uarousal, states names roughnrem.
    %         Transitions are found on the fine bouts and controls searched inside
    %         the coarse ones

    properties (Dependent)
        focus   % 1 x N cell  the transition labels, column 3 of param.merge. The
                % subset of eventlist a PIV run walks; derived, never set
    end
    %   info, field for field. Nothing here is an echo of an argument.
    %     dropped    struct   empty_band / out_of_range, per event_pick_excursions
    %     clipped    1 x K    events whose span hit the trace ends
    %     sign_ok    bool     every event moved the way its polarity prescribed
    %     sign_bad   1 x K    the ones that did not
    %     P          struct   what the picking derived: sg_win_s, peaktol, fps
    %     n_pair     int      transition pairs merge_bouts produced
    %     quiet      struct   counts per state, n_bout, and the matched-pair cap

    methods
        function obj = analysis_event(diameter, t_axis)
        % IN   diameter    1 x T float   caliber series in PX. An ABSOLUTE field:
        %                                the *changes fields are median-subtracted,
        %                                and the amplitude banding needs a caliber
        %      t_axis      1 x T float   frame times (s)
            arguments
                diameter (1,:) double {mustBeNonempty}
                t_axis   (1,:) double {mustBeNonempty}
            end
            if numel(diameter) ~= numel(t_axis)
                error('analysis_event:lengthMismatch', ...
                    'diameter has %d samples, t_axis has %d.', ...
                    numel(diameter), numel(t_axis));
            end
            obj.diameter = diameter;
            obj.t_axis   = t_axis;
        end

        function v = get.focus(obj)
            if isempty(obj.param.merge)
                v = {};
                return
            end
            v = obj.param.merge(:, 3)';
        end

        function arousaltransition(obj, state_bout, merge, opt)
        % AROUSALTRANSITION  Find the endpoints of every state-transition excursion.
        %   Fills smooth_diameter, eventlist and info.
        %
        % IN   state_bout struct       [nBout x 2] frame-index tables per state, from
        %                              state_image.state_idx. This class knows nothing
        %                              about sleep scoring, so they come in
        %      merge      N x 4 cell   {earlierState, laterState, label, polarity}.
        %                              Omitted = obj.param.merge; given, it is
        %                              WRITTEN there, so the object and the run
        %                              cannot disagree about which partition ran
        %      search_s   float s      search window each side of a transition
        %      sg_win_s   float s      Savitzky-Golay window of the picking trace
        %      peaktol    float 0..1   band level, median toward the extreme
        %
        %   The two primitives are called here in order rather than behind one
        %   helper: merge_bouts reads the bout tables and never the trace,
        %   pick_excursions reads the trace and never the bouts. Which one sees what
        %   is why the polarity is prescribed rather than measured.
            arguments
                obj
                state_bout     struct
                merge          cell = {}
                opt.search_s   (1,1) double = 30
                opt.sg_win_s   (1,1) double = 3
                opt.peaktol    (1,1) double = 0.80
            end
            if ~isempty(merge)
                obj.param.merge = merge;
            end
            if isempty(obj.param.merge)
                error('analysis_event:noMerge', ...
                    ['no transition table. Pass one, or set param.merge -- ' ...
                     '{earlier, later, label, polarity}, one row per transition.']);
            end

            % 1. Structural: which stretches to search and what polarity to expect
            pairs = event_merge_bouts(state_bout, obj.param.merge, obj.t_axis, opt.search_s);

            % 2. Where on the trace each endpoint actually lands
            [obj.eventlist, pickinfo] = event_pick_excursions( ...
                obj.diameter, obj.t_axis, pairs, opt.sg_win_s, opt.peaktol);

            obj.smooth_diameter = pickinfo.dsg;
            pickinfo            = rmfield(pickinfo, 'dsg');
            pickinfo.n_pair     = numel(pairs);
            obj.info            = pickinfo;

            % 3. One line where it happens; info holds the rest
            drop = obj.info.dropped;
            fprintf(['%d pairs -> %d events (px) | dropped %d empty band, ' ...
                     '%d out of range | %d clipped\n'], ...
                numel(pairs), numel(obj.eventlist), ...
                drop.empty_band, drop.out_of_range, numel(obj.info.clipped));
            if ~obj.info.sign_ok
                warning('analysis_event:signMismatch', ...
                    'event(s) %s contradict their prescribed polarity', ...
                    mat2str(obj.info.sign_bad));
            end
        end

        function control(obj, state_bout, opt)
        % CONTROL  Add the stable stretches to eventlist, as (state, none) rows.
        %   A control is only worth having next to the event it controls, so these
        %   go into the same list and everything downstream walks both.
        %
        % IN   state_bout  struct       the SAME tables arousaltransition was given.
        %                               It found the boundaries BETWEEN states; this
        %                               searches INSIDE them, so a different
        %                               partition would make the avoidance meaningless
        %      states      cell         which state_bout fields to search, {} = all
        %      search_frac 1 x 2 float  the portion of each bout to scan
        %      dd_frac     float        how far the two ENDS may differ, as a
        %                               fraction of the median |diameter_change| of
        %                               the transitions this object already found.
        %                               0.1 = the pair must match to a tenth of what
        %                               an event moves
        %      range_frac  float        optional cap on the excursion BETWEEN the
        %                               ends, same scaling. Inf = report, do not gate
        %      len_frac    1 x K float  window length per bout, see event_findquiet
        %      len_min_s   float s      shortest window worth reporting
        %      margin      int frames   padding around each transition a window may
        %                               not touch. The PIV object cuts
        %                               from-halfwin : to+halfwin, so a control
        %                               inside that margin shares frames with the
        %                               event it is a control for
        %      seed        int          which windows the random pick returns
        %
        %   dd_frac is a fraction and not an absolute number so the cap carries no
        %   unit of its own and says what it means.
        %
        %   range and sd are backfilled onto the TRANSITION rows too, over their own
        %   from:to. One struct array, one field set, one meaning per field.
            arguments
                obj
                state_bout      struct
                opt.states      cell = {}
                                     % omitted = obj.param.states; given, written there
                opt.search_frac (1,2) double = [1/3 2/3]
                opt.dd_frac     (1,1) double {mustBePositive} = 0.1
                opt.range_frac  (1,1) double {mustBePositive} = Inf
                opt.len_frac    (1,:) double = [1 0.5 0.25]
                opt.len_min_s   (1,1) double = 2
                opt.margin      (1,1) double = 2
                opt.seed        (1,1) double = 0
            end
            if ~isempty(opt.states)
                obj.param.states = opt.states;
            end
            if isempty(obj.eventlist) || isempty(obj.smooth_diameter)
                error('analysis_event:noTransitions', ...
                    ['run arousaltransition first; controls share its smoothed ' ...
                     'trace and avoid its events.']);
            end

            % 1. Back to the transitions only, so a second call does not stack
            keep     = ~strcmp({obj.eventlist.pol}, 'none');
            eventrow = obj.eventlist(keep);
            dsg      = obj.smooth_diameter;

            % 2. Frames a stable window may not touch: every transition plus the
            %    margin the PIV object will cut around it
            n_frame  = numel(obj.diameter);
            excluded = false(1, n_frame);
            for e = 1:numel(eventrow)
                span_from = max(1, eventrow(e).from - opt.margin);
                span_to   = min(n_frame, eventrow(e).to + opt.margin);
                excluded(span_from:span_to) = true;
            end

            % 3. How closely the two ends must match, set by the events themselves
            dD_med    = median(abs([eventrow.diameter_change]));
            max_dd    = opt.dd_frac    * dD_med;
            max_range = opt.range_frac * dD_med;

            % 4. The stable stretches
            [quietrow, quietinfo] = event_findquiet(obj.diameter, obj.t_axis, ...
                state_bout, states = obj.param.states, ...
                search_frac = opt.search_frac, max_dd = max_dd, ...
                max_range = max_range, len_frac = opt.len_frac, ...
                len_min_s = opt.len_min_s, exclude = excluded, dsg = dsg, ...
                seed = opt.seed);

            % 5. Same fields on both, then one chronological list
            for e = 1:numel(eventrow)
                seg = eventrow(e).from : eventrow(e).to;
                eventrow(e).range = max(dsg(seg)) - min(dsg(seg));
                eventrow(e).sd    = std(dsg(seg));
            end
            if isempty(quietrow)
                obj.eventlist = eventrow;
            else
                merged   = [eventrow, orderfields(quietrow, eventrow)];
                [~, ord] = sort([merged.from]);
                obj.eventlist = merged(ord);
            end
            obj.info.quiet = quietinfo;

            % 6. One line where it happens. n is BOUTS, not windows -- the lengths
            %    from one bout share its tissue even when they no longer share frames
            fprintf(['%d transitions + %d matched pairs from %d bouts | ' ...
                     'ends must agree to %.3f px (%.0f%% of the median event)\n'], ...
                numel(eventrow), numel(quietrow), quietinfo.n_bout, max_dd, ...
                100*opt.dd_frac);
            for s = string(fieldnames(quietinfo.counts))'
                c = quietinfo.counts.(s);
                fprintf(['  %-12s %2d bouts -> %2d windows | %d bout too short, ' ...
                         '%d length had no pair\n'], ...
                    s, c.bouts, c.windows, c.too_short, c.no_matched_pair);
            end
        end

        function findexcursions(obj, state_bout, opt)
        % FINDEXCURSIONS  Every diameter excursion inside every state bout, once.
        %   Fills obj.excursion. It SELECTS NOTHING -- a row is any pair of
        %   consecutive turning points, and the properties a selection would key
        %   on are columns. A rule that would have been a re-run is a filter:
        %
        %     T(T.state=="nrem" & T.pol=="dilation" & T.med2max_px > 0, :)
        %     T(T.prom_min_px >= 2 & T.dur_s >= 10, :)
        %     T(isnan(T.contained_by), :)
        %
        %   so two results that disagree differ by a filter that is written down.
        %   see design/EXCURSION_DESIGN.md for why, and CLAUDE_LOG.md for the day
        %   that produced five incomparable event sets from five re-runs.
        %
        % IN   state_bout   struct       [nBout x 2] frame-index bout tables per state
        %      states       cell         which state_bout fields to walk. {} = all
        %      prom_px      1 x P float  MinPeakProminence values to sweep, in the
        %                   unit of diameter. A long rise split by a ripple is one
        %                   row at high prominence and several at low, and both are
        %                   kept -- prom_min_px says which
        %      edge_margin  int frames   clearance an endpoint keeps from the bout
        %                   edge. Required: it must cover whatever half-window a
        %                   later stage puts around the endpoints, and this class
        %                   cannot know that number
        %      min_frame    int          shortest span kept, in frames
        %      upper_edge   1 x T float  the vessel's two wall positions, for
        %      lower_edge                edge_r. [] = leave edge_r NaN
        % OUT  obj.excursion, height x 17 table
        %
        %   MED2MIN AND MED2MAX ARE THE CLASSIFICATION. They are the two extremes
        %   of the excursion measured from the bout's median calibre, so a rise
        %   from 3 px below rest to 5 px above is (-3, +5) and a fall over the same
        %   ground is also (-3, +5) with pol reversed. Direction lives in pol and
        %   geometry lives in these two, and neither contaminates the other. There
        %   is no category column: med2min >= 0 is an excursion that stayed above
        %   rest, med2max <= 0 one that stayed below, and anything else crosses.
            arguments
                obj
                state_bout       struct
                opt.states       cell = {}
                opt.prom_px      (1,:) double {mustBeNonnegative} = [0 1 2 3]
                opt.edge_margin  (1,1) double {mustBeInteger, mustBeNonnegative}
                opt.min_frame    (1,1) double {mustBeInteger, mustBePositive} = 3
                opt.upper_edge   (1,:) double = []
                opt.lower_edge   (1,:) double = []
                opt.verbose      (1,1) logical = true
            end
            if isempty(obj.smooth_diameter)
                error('analysis_event:noSmooth', ...
                    ['smooth_diameter is empty; run arousaltransition first. The ' ...
                     'excursions have to be picked off the SAME trace the ' ...
                     'transitions were, or the two are not comparable.']);
            end
            n_sample = numel(obj.diameter);
            states = opt.states;
            if isempty(states)
                states = fieldnames(state_bout)';
            end

            % 1. Walk every bout at every prominence. Duplicates are expected and
            %    are resolved afterwards, because the same rise found at two
            %    prominences is one excursion with two witnesses
            raw = zeros(0, 8);   % from to state_i base med2min med2max prom bout
            for s = 1:numel(states)
                name = states{s};
                if ~isfield(state_bout, name)
                    continue
                end
                bouts = state_bout.(name);
                for i = 1:size(bouts, 1)
                    b1 = max(1, bouts(i,1));
                    b2 = min(n_sample, bouts(i,2));
                    s1 = b1 + opt.edge_margin;
                    s2 = b2 - opt.edge_margin;
                    if s2 - s1 < opt.min_frame
                        continue
                    end
                    base = median(obj.diameter(b1:b2));
                    seg  = obj.smooth_diameter(s1:s2);
                    for p = opt.prom_px
                        [loc, ~] = local_turningpoints(seg, p);
                        loc = loc + s1 - 1;
                        for k = 1:numel(loc) - 1
                            f = loc(k);
                            t = loc(k+1);
                            if t - f < opt.min_frame
                                continue
                            end
                            span = f:t;
                            raw(end+1,:) = [f, t, s, base, ...
                                min(obj.diameter(span)) - base, ...
                                max(obj.diameter(span)) - base, p, i]; %#ok<AGROW>
                        end
                    end
                end
            end
            if isempty(raw)
                obj.excursion = table();
                return
            end

            % 2. One row per (from, to), carrying BOTH ends of the prominence range
            %    it was found over. The two say different things and one alone
            %    conflates two kinds of excursion:
            %      prom_min 0, prom_max 0      a ripple; asking for 1 removes it
            %      prom_min 0, prom_max high   a clean rise with no ripple inside it,
            %                                  so the same two turning points at
            %                                  every scale
            %      prom_min high               a coarse rise that EXISTS only once
            %                                  the ripples inside it are ignored
            [~, ord] = sortrows(raw, [1 2 7]);
            raw = raw(ord, :);
            [pairs, keep] = unique(raw(:, [1 2]), 'rows', 'first');
            prom_max = zeros(size(pairs, 1), 1);
            for a = 1:size(pairs, 1)
                same = raw(:,1) == pairs(a,1) & raw(:,2) == pairs(a,2);
                prom_max(a) = max(raw(same, 7));
            end
            [keep, back] = sort(keep);
            prom_max = prom_max(back);
            raw = raw(keep, :);
            n = size(raw, 1);

            % 3. Containment and overlap. Without these the table cannot say what
            %    a count of its rows is a count OF
            contained_by = nan(n, 1);
            overlaps_n   = zeros(n, 1);
            for a = 1:n
                for b = 1:n
                    if a == b
                        continue
                    end
                    inside = raw(b,1) <= raw(a,1) && raw(b,2) >= raw(a,2);
                    if inside
                        wider = (raw(b,2) - raw(b,1)) > (raw(a,2) - raw(a,1));
                        if wider && (isnan(contained_by(a)) || ...
                                (raw(b,2)-raw(b,1)) < (raw(contained_by(a),2)-raw(contained_by(a),1)))
                            contained_by(a) = b;      % the TIGHTEST container
                        end
                    end
                    if ~(raw(b,2) <= raw(a,1) || raw(b,1) >= raw(a,2))
                        overlaps_n(a) = overlaps_n(a) + 1;
                    end
                end
            end

            % 4. Per-row quantities the walk did not need
            from   = raw(:,1);
            to     = raw(:,2);
            state  = string(states(raw(:,3))');
            pol    = repmat("constriction", n, 1);
            dD_px  = obj.smooth_diameter(to)' - obj.smooth_diameter(from)';
            pol(dD_px > 0) = "dilation";
            dD_raw_px = obj.diameter(to)' - obj.diameter(from)';
            dur_s  = obj.t_axis(to)' - obj.t_axis(from)';
            range_px = nan(n,1);
            sd_px    = nan(n,1);
            edge_r   = nan(n,1);
            has_edge = ~isempty(opt.upper_edge) && ~isempty(opt.lower_edge);
            for a = 1:n
                span = from(a):to(a);
                range_px(a) = max(obj.smooth_diameter(span)) - min(obj.smooth_diameter(span));
                sd_px(a)    = std(obj.smooth_diameter(span));
                if has_edge
                    u = opt.upper_edge(span);
                    v = opt.lower_edge(span);
                    if std(u) > 0 && std(v) > 0
                        edge_r(a) = corr(u(:), v(:));
                    end
                end
            end

            obj.excursion = table(from, to, dur_s, state, pol, ...
                raw(:,8), raw(:,4), raw(:,5), raw(:,6), ...
                dD_px, dD_raw_px, range_px, sd_px, ...
                raw(:,7), prom_max, contained_by, overlaps_n, edge_r, ...
                'VariableNames', {'from','to','dur_s','state','pol', ...
                    'bout_idx','base_px','med2min_px','med2max_px', ...
                    'dD_px','dD_raw_px','range_px','sd_px', ...
                    'prom_min_px','prom_max_px','contained_by','overlaps_n','edge_r'});
            obj.excursion.Properties.VariableUnits = { ...
                'frame','frame','s','','','','px','px','px', ...
                'px','px','px','px','px','px','','',''};
            obj.excursion = sortrows(obj.excursion, 'from');
            % contained_by held indices into the pre-sort order; carry them over
            obj.excursion.contained_by = repoint(contained_by, from, obj.excursion.from);

            obj.info.excursion = struct('states', {states}, 'prom_px', opt.prom_px, ...
                'edge_margin', opt.edge_margin, 'min_frame', opt.min_frame, ...
                'n_maximal', nnz(isnan(contained_by)));
            if opt.verbose
                fprintf(['excursion  %d rows from %d bouts | %d maximal, ' ...
                         '%d contained | prominence %s\n'], n, ...
                    size(unique(raw(:,[3 8]), 'rows'), 1), nnz(isnan(contained_by)), ...
                    nnz(~isnan(contained_by)), mat2str(opt.prom_px));
            end
        end

        function save2disk(obj, name, savepath)
            analysis_event = obj;
            save(fullfile(savepath, [name, '.mat']), 'analysis_event')
        end
    end
end

% ---------------------------------------------------------------------------
function [loc, is_max] = local_turningpoints(y, prom)
%LOCAL_TURNINGPOINTS  Maxima and minima of y, interleaved, in y's own indices.
%   prom = 0 asks findpeaks for every local extremum; above 0 it ignores a bump
%   that never gets that far from its neighbours, which is what stops a ripple
%   partway up a long rise from ending it.
    if prom > 0
        [~, lmax] = findpeaks(y,  'MinPeakProminence', prom);
        [~, lmin] = findpeaks(-y, 'MinPeakProminence', prom);
    else
        [~, lmax] = findpeaks(y);
        [~, lmin] = findpeaks(-y);
    end
    loc    = [lmax(:).', lmin(:).'];
    is_max = [true(1, numel(lmax)), false(1, numel(lmin))];
    [loc, ord] = sort(loc);
    is_max = is_max(ord);
end

% ---------------------------------------------------------------------------
function out = repoint(idx, from_old, from_new)
%REPOINT  Carry a column of row indices through a re-sort of the rows.
    out = nan(size(idx));
    have = ~isnan(idx);
    target = from_old(idx(have));
    [~, out(have)] = ismember(target, from_new);
    out(out == 0) = NaN;
end
