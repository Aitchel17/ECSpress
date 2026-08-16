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

        function save2disk(obj, name, savepath)
            analysis_event = obj;
            save(fullfile(savepath, [name, '.mat']), 'analysis_event')
        end
    end
end
