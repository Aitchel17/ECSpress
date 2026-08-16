classdef vfield_profile < handle
%VFIELD_PROFILE  Radial profiles from many fields, compared against quiet rows.
%   One row per event. Each row holds the n_bin x n_wedge block that
%   vfield_polar.accumulate produced, plus what it needs to be believed: polarity,
%   state, whether the field was gated, and how far each wedge reached. The
%   displacement fields themselves are not kept, so a whole cohort fits in memory.
%
%   THIS CLASS EXISTS BECAUSE OF ONE BUG.
%   Averaging over a set whose membership changes with the coordinate. It happened
%   three times in one day -- a ring mean over wedges that differed per ring, a
%   table whose sector set changed per row, and a radial rise that was half
%   composition. The two rules that came out of it are enforced here as method
%   signatures rather than left to the caller:
%
%     1  a wedge is in or out for the WHOLE comparison interval, never per bin
%     2  only cells where every row in the comparison has data, and no omitnan
%        stepping over the ones where it does not
%
%   gate_support() computes the verdict, and curve() and difference() will not
%   run without a bin range to compute it over.
%
%   THE CONTROL IS event_findquiet's OWN ROWS.
%   Not a window built for the occasion. The quiet rows share the geometry, the
%   optics and the recording, so passing them through this identical pipeline is
%   the only thing that says whether a number is real. Rows with pol == "none".
%
%   POLARITY IS NOT FOLDED BY DEFAULT.
%   Displacement reverses with polarity because that is geometry, but a VOLUME
%   change need not. Folding cancels anything that comes out the same sign both
%   ways -- which is rectification, one of the things this is here to look for.
%   See CLAUDE_LOG.md
%
%   see also VFIELD_POLAR, EVENT_FINDQUIET

    properties
        param = struct( ...
            'min_row_cell', 3, ...     % int   rows a cell needs before it counts
            'verbose',      true)      % bool  one line per append
    end

    properties (SetAccess = private)
        % settled at construction : every row is expressed on this axis, so a row
        % measured on another one cannot be appended
        bin_edges_um   % 1 x nB+1 float  radial bin edges, MICRONS from the wall
        n_wedge        % int             angular bins

        rows = struct([])  % 1 x n_row struct. DATA, one row per event:
                           %   event_idx  int
                           %   pol        string   "dilation"/"constriction"/"none"
                           %   state      string   the transition
                           %   dD         float um diameter change
                           %   sgn        int      +1 / -1
                           %   gated      bool     provenance, carried from the field
                           %   cells      struct   n_bin x n_wedge, see append
                           %   reach_bin  1 x n_wedge float
                           %   mag_kept   n_bin x 1 float um   [] if not audited
                           %   mag_cut    n_bin x 1 float um   [] if not audited
    end

    properties (Dependent)
        bin_center_um  % 1 x nB float
        n_bin          % int
        n_row          % int
    end

    methods
        function obj = vfield_profile(bin_edges_um, n_wedge)
            arguments
                bin_edges_um (1,:) double
                n_wedge (1,1) double {mustBePositive, mustBeInteger}
            end
            obj.bin_edges_um = bin_edges_um;
            obj.n_wedge = n_wedge;
        end

        function v = get.bin_center_um(obj)
            lo = obj.bin_edges_um(1:end-1);
            hi = obj.bin_edges_um(2:end);
            v = (lo + hi) / 2;
        end

        function v = get.n_bin(obj)
            v = numel(obj.bin_edges_um) - 1;
        end

        function v = get.n_row(obj)
            v = numel(obj.rows);
        end

        function append(obj, cells, opt)
        %APPEND  Add one event's block. Takes the BLOCK, not a vfield_polar.
        %   Deliberately not the object: if this took the field it would have to
        %   drive its stages, and the order of measure -> gate_wedge -> accumulate
        %   belongs to vfield_polar. The caller runs them and hands over what came
        %   out, so a half-run field cannot be silently finished here.
        % IN  cells  struct of n_bin x n_wedge, from vfield_polar.accumulate
            arguments
                obj
                cells struct
                opt.event_idx (1,1) double
                opt.pol (1,1) string
                opt.state (1,1) string
                opt.dD (1,1) double
                opt.gated (1,1) logical
                opt.mag_kept double = []
                opt.mag_cut double = []
            end
            got = size(cells.divergence);
            if ~isequal(got, [obj.n_bin, obj.n_wedge])
                error('vfield_profile:axisMismatch', ...
                    'block is %d x %d, this profile is %d x %d.', ...
                    got(1), got(2), obj.n_bin, obj.n_wedge);
            end
            if obj.n_row > 0 && obj.rows(1).gated ~= opt.gated
                error('vfield_profile:gateMismatch', ...
                    'row %d is gated=%d, this profile holds gated=%d.', ...
                    opt.event_idx, opt.gated, obj.rows(1).gated);
            end

            sgn = sign(opt.dD);
            if sgn == 0
                sgn = 1;
            end
            row = struct('event_idx', opt.event_idx, 'pol', opt.pol, ...
                'state', opt.state, 'dD', opt.dD, 'sgn', sgn, ...
                'gated', opt.gated, 'cells', cells, ...
                'reach_bin', cells.reach_bin, ...
                'mag_kept', opt.mag_kept, 'mag_cut', opt.mag_cut);
            if obj.n_row == 0
                obj.rows = row;
            else
                obj.rows(end+1) = row;
            end
            if obj.param.verbose
                fprintf('append        ev %d  %s  %s  dD %+.2f\n', ...
                    opt.event_idx, opt.pol, opt.state, opt.dD);
            end
        end

        function v = peak(obj, quantity, rowsel)
        %PEAK  The largest magnitude one quantity reaches anywhere in the cohort.
        %   A cohort number, which is why it is here: max over ROWS is the same
        %   shape of question as wedges_over and widest_common, and a single row
        %   cannot answer it.
        % IN  quantity  char, which field of cells to look at
        %     rowsel    1 x n_row bool, default every row
        % OUT v         1 x 1 float, max |value| over the selected rows. 1 when
        %               nothing finite was found, so a caller can divide by it
        %
        %   why  the display layer wants ONE limit across every panel -- pass a
        %        per-event limit and each map self-normalises, so a weak event
        %        and a strong one come out looking the same
        %        (showpiv.plot_overlay says the same thing from its side)
        %   note named for what it measures, not for what it is for. This class
        %        does not know that the number becomes a colour axis, and it
        %        should not -- see CLAUDE.md, engines stay pure
            arguments
                obj
                quantity (1,:) char
                rowsel   (1,:) logical = true(1, obj.n_row)
            end
            v = 0;
            for k = find(rowsel)
                M = obj.rows(k).cells.(quantity);
                m = max(abs(M(:)), [], 'omitnan');
                if isfinite(m)
                    v = max(v, m);
                end
            end
            if v <= 0
                v = 1;
            end
        end

        function support = gate_support(obj, rowsel, quantity)
        %GATE_SUPPORT  Which cells every row in the comparison actually has.
        %   n_bin x n_wedge, NOT n_bin x n_wedge x n_row. There is one partition,
        %   so there is one common support; a per-row answer would only invite
        %   someone to average over a set that changes.
        % IN  rowsel    1 x n_row bool, the rows that must ALL have data
        %     quantity  char, which field of cells to test
            arguments
                obj
                rowsel (1,:) logical
                quantity (1,:) char = 'divergence'
            end
            support = true(obj.n_bin, obj.n_wedge);
            idx = find(rowsel);
            for k = idx
                M = obj.rows(k).cells.(quantity);
                support = support & isfinite(M);
            end
        end

        function [wedge_ok, support] = wedges_over(obj, rowsel, bin_range, quantity)
        %WEDGES_OVER  Wedges with data at every bin of the interval, for every row.
        %   Rule 1. A wedge is in or out for the WHOLE interval; there is no
        %   per-bin membership, because that is what made the ring mean wrong.
        % IN  bin_range 1 x 2 int, [first last] bin index inclusive
        % OUT wedge_ok  1 x n_wedge bool
            arguments
                obj
                rowsel (1,:) logical
                bin_range (1,2) double
                quantity (1,:) char = 'divergence'
            end
            support = obj.gate_support(rowsel, quantity);
            span = bin_range(1):bin_range(2);
            wedge_ok = all(support(span, :), 1);
        end

        function [bin_range, wedge_ok] = widest_common(obj, rowsel, quantity)
        %WIDEST_COMMON  The longest bin span some wedge covers for every row.
        %   Every comparison needs an interval, and picking one by hand at the
        %   call site is how a magic number gets into a figure. This reports the
        %   longest one the data actually supports; a caller wanting a narrower
        %   one still passes it explicitly.
        %
        %   IT MAXIMISES THE SPAN, NOT THE WEDGE COUNT, so it is the right default
        %   for a PLOT and the wrong one for a CLAIM. Asking a wedge to reach the
        %   outer edge throws away wedges that stop sooner, and a claim about the
        %   inner radii does not need them to. For a statement at radius r, pass
        %   [1 bin(r)] and keep the wedges that reach only that far.
        % OUT bin_range 1 x 2 int, [first last]. [NaN NaN] if nothing is common
            arguments
                obj
                rowsel (1,:) logical
                quantity (1,:) char = 'divergence'
            end
            support = obj.gate_support(rowsel, quantity);
            bin_range = [NaN NaN];
            wedge_ok = false(1, obj.n_wedge);
            best = 0;
            for first = 1:obj.n_bin
                for last = first:obj.n_bin
                    ok = all(support(first:last, :), 1);
                    if ~any(ok)
                        break
                    end
                    span = last - first + 1;
                    if span > best
                        best = span;
                        bin_range = [first last];
                        wedge_ok = ok;
                    end
                end
            end
        end

        function [curves, wedge_ok] = curve_rows(obj, rowsel, bin_range, opt)
        %CURVE_ROWS  One curve per row, over a fixed wedge set. No omitnan.
        % OUT curves    n_selected x n_bin float, NaN outside the interval
        %     wedge_ok  1 x n_wedge bool, the set every curve used
            arguments
                obj
                rowsel (1,:) logical
                bin_range (1,2) double
                opt.quantity (1,:) char = 'volume_out'
                opt.fold (1,1) logical = false
            end
            wedge_ok = obj.wedges_over(rowsel, bin_range, opt.quantity);
            idx = find(rowsel);
            curves = nan(numel(idx), obj.n_bin);
            span = bin_range(1):bin_range(2);
            for i = 1:numel(idx)
                k = idx(i);
                M = obj.rows(k).cells.(opt.quantity);
                block = M(span, wedge_ok);
                % plain mean, not omitnan : wedges_over already guaranteed every
                % one of these is finite, so a NaN here is a fault to see
                curves(i, span) = mean(block, 2);
                if opt.fold
                    curves(i, span) = curves(i, span) * obj.rows(k).sgn;
                end
            end
        end

        function out = difference(obj, rowsel, bin_range, opt)
        %DIFFERENCE  A group against the quiet rows, through the identical path.
        %   The quiet rows enter wedges_over too, so both arms use one wedge set.
        % OUT out  table : r_um, group, quiet, diff, sigma, z
            arguments
                obj
                rowsel (1,:) logical
                bin_range (1,2) double
                opt.quantity (1,:) char = 'volume_out'
                opt.fold (1,1) logical = false
            end
            isq = obj.quiet_rows();
            % A difference against nothing is not a difference. Without this the
            % table comes back all NaN and only the verbose line says why, which is
            % how a whole run of transitions-only ended up looking like a result
            if ~any(isq)
                error('vfield_profile:noQuietRows', ...
                    ['no pol=="none" rows in this profile. The controls have to go ' ...
                     'through PIV too -- see piv_integration_testbed cell 11.']);
            end
            both = rowsel | isq;
            wedge_ok = obj.wedges_over(both, bin_range, opt.quantity);
            if ~any(wedge_ok)
                error('vfield_profile:noCommonSupport', ...
                    'no wedge has data across bins %d-%d for every row.', ...
                    bin_range(1), bin_range(2));
            end

            g = obj.group_curve(rowsel, wedge_ok, bin_range, opt.quantity, opt.fold);
            q = obj.group_curve(isq, wedge_ok, bin_range, opt.quantity, false);
            diff_val = g.mean - q.mean;
            sigma = sqrt(g.sem.^2 + q.sem.^2);
            z = diff_val ./ sigma;

            span = bin_range(1):bin_range(2);
            out = table(obj.bin_center_um(span)', g.mean(span)', q.mean(span)', ...
                diff_val(span)', sigma(span)', z(span)', ...
                'VariableNames', {'r_um','group','quiet','diff','sigma','z'});
            if obj.param.verbose
                fprintf('difference    %d rows vs %d quiet, %d/%d wedges\n', ...
                    nnz(rowsel), nnz(isq), nnz(wedge_ok), obj.n_wedge);
            end
        end

        function out = increment(obj, rowsel, bin_range, opt)
        %INCREMENT  Bin-to-bin change of the quiet-subtracted curve.
        %   A cumulative that has stopped rising is not the same claim as one that
        %   is still rising, and the level cannot tell them apart.
            arguments
                obj
                rowsel (1,:) logical
                bin_range (1,2) double
                opt.quantity (1,:) char = 'volume_out'
                opt.fold (1,1) logical = false
            end
            d = obj.difference(rowsel, bin_range, quantity = opt.quantity, fold = opt.fold);
            step = diff(d.diff);
            % the two levels share their quiet arm, so their errors are correlated;
            % this is the conservative bound, not the exact one
            step_sigma = sqrt(d.sigma(1:end-1).^2 + d.sigma(2:end).^2);
            out = table(d.r_um(2:end), step, step_sigma, step ./ step_sigma, ...
                'VariableNames', {'r_um','step','sigma','z'});
        end

        function f = frac_absorbed(obj, rowsel, bin_range, opt)
        %FRAC_ABSORBED  1 - volume_out(outer) / volume_out(inner). A method, not a
        %   property: it is a question about two radii the caller picks, and
        %   storing it would hide which two.
            arguments
                obj
                rowsel (1,:) logical
                bin_range (1,2) double
                opt.fold (1,1) logical = false
            end
            d = obj.difference(rowsel, bin_range, quantity = 'volume_out', fold = opt.fold);
            inner = d.diff(1);
            outer = d.diff(end);
            f = 1 - outer / inner;
        end

        function tbl = audit_gate(obj, rowsel)
        %AUDIT_GATE  Pool the per-row gate audits. Empty if none were supplied.
        %   mag_cut > mag_kept means the gates took the LARGER displacements at
        %   that radius, which no quiet subtraction can undo: the bias scales with
        %   amplitude and the quiet rows have none. See CLAUDE_LOG.md
            arguments
                obj
                rowsel (1,:) logical
            end
            idx = find(rowsel);
            keep = nan(obj.n_bin, numel(idx));
            cut = nan(obj.n_bin, numel(idx));
            for i = 1:numel(idx)
                k = idx(i);
                if isempty(obj.rows(k).mag_kept)
                    continue
                end
                keep(:, i) = obj.rows(k).mag_kept;
                cut(:, i) = obj.rows(k).mag_cut;
            end
            mk = median(keep, 2, 'omitnan');
            mc = median(cut, 2, 'omitnan');
            tbl = table(obj.bin_center_um(:), mk, mc, mc ./ mk, ...
                'VariableNames', {'r_um','mag_kept_um','mag_cut_um','cut_over_kept'});
        end

        function sel = quiet_rows(obj)
        %QUIET_ROWS  event_findquiet's own rows, as they arrive in events.
            sel = false(1, obj.n_row);
            for k = 1:obj.n_row
                sel(k) = obj.rows(k).pol == "none";
            end
        end

        function sel = select(obj, opt)
        %SELECT  A row mask by polarity, state or amplitude. Never includes quiet.
            arguments
                obj
                opt.pol string = string.empty
                opt.state string = string.empty
                opt.min_abs_dD (1,1) double = 0
            end
            sel = ~obj.quiet_rows();
            for k = 1:obj.n_row
                if ~sel(k)
                    continue
                end
                if ~isempty(opt.pol) && ~any(obj.rows(k).pol == opt.pol)
                    sel(k) = false;
                end
                if ~isempty(opt.state) && ~any(obj.rows(k).state == opt.state)
                    sel(k) = false;
                end
                if abs(obj.rows(k).dD) < opt.min_abs_dD
                    sel(k) = false;
                end
            end
        end

        function s = to_struct(obj)
        %TO_STRUCT  Plain struct for save(). A class in a .mat needs the classdef.
            s = struct('bin_edges_um', obj.bin_edges_um, 'n_wedge', obj.n_wedge, ...
                'param', obj.param, 'rows', obj.rows);
        end
    end

    methods (Static)
        function obj = from_struct(s)
        %FROM_STRUCT  Rebuild from to_struct output.
            obj = vfield_profile(s.bin_edges_um, s.n_wedge);
            obj.param = s.param;
            obj.rows = s.rows;
        end
    end

    methods (Access = private)
        function out = group_curve(obj, rowsel, wedge_ok, bin_range, quantity, fold)
        % Mean and standard error across rows, on a wedge set the caller fixed.
        % Every row uses the same wedges and the same bins, so what varies between
        % rows is the measurement and nothing else
            idx = find(rowsel);
            span = bin_range(1):bin_range(2);
            per_row = nan(numel(idx), obj.n_bin);
            for i = 1:numel(idx)
                k = idx(i);
                M = obj.rows(k).cells.(quantity);
                block = M(span, wedge_ok);
                per_row(i, span) = mean(block, 2);
                if fold
                    per_row(i, span) = per_row(i, span) * obj.rows(k).sgn;
                end
            end
            out = struct();
            out.mean = mean(per_row, 1, 'omitnan');
            out.sem = std(per_row, 0, 1, 'omitnan') ./ sqrt(sum(isfinite(per_row), 1));
            out.n = sum(isfinite(per_row), 1);
        end
    end
end
