classdef piv_figure < handle
%PIV_FIGURE  Display layer for per-event ensemble-PIV results.
%   Holds the two settings that must stay constant across events:
%     quiver scale   arrow length = displacement(px) * frame_span * scale, and
%                    frame_span is pinned to 1 because uv is already a total
%                    from->to displacement (see piv_run_events)
%     colour limits  computed once over all events, symmetric about 0, so the maps
%                    show amplitude differences instead of self-normalising
%
%   THE MAPS ARE PAINTED, NOT MEASURED.
%   vfield_polar does not emit a dense map; it bins triangles. paint() spreads a
%   per-cell value back over the frame so it can be looked at where it came from.
%   Nothing downstream reads a painted map.
%
%   radprofile() USED TO BE WRONG.
%   It ran mean(cells, 2, 'omitnan') -- a mean over the wedges present at each
%   ring, and a wedge that runs off the frame stops sooner than its neighbours,
%   so the inner radii averaged more wedges than the outer ones. It now goes
%   through vfield_profile, which fixes one wedge set across the whole interval.
%   See CLAUDE_LOG.md
%
% IN   run      1xN struct from piv_run_events
%               (id,state,pol,from,to,rise_s,diameter_change,nfr,uv,corr)
%      fields   1xN cell of vfield_polar, from piv_polar_events
%      profile  vfield_profile from the same call, one row per event
%      ctx      struct: S, coremask, exclmask, t_axis, dtrace, dsg, fps, p2u,
%               halfwin, events, and optionally state_idx / states
% OUT  obj      handle object; call the methods below
%
%   USAGE
%     [P, V] = piv_polar_events(piv_run, coremask, p2u, gated = false);
%     F = piv_figure(piv_run, V, P, ctx);
%     F.summary();                     % text table
%     F.map(2, 'divergence');          % overlay + quiver for event index 2
%     F.panel(2);                      % four panels, one figure
%     F.tiles('dilation');
%     F.radprofile();  F.signcheck();  F.slices(2);

    properties
        run                     % 1 x N struct   per-event PIV results
        fields                  % 1 x N cell     vfield_polar, one per event
        profile                 % vfield_profile the cohort, one row per event
        ctx                     % struct         stack, masks, traces, geometry
        lims                    % struct         symmetric limits shared by all events
        scale   = 5             % float          quiver length multiplier
        alpha   = 0.55          % float          overlay transparency
        cmap    = 'turbo'       % char
        headsz  = 5             % float px       absolute arrowhead size
        linew   = 1             % float
    end

    properties (Constant, Access = private)
        req_fields = {'S','coremask','exclmask','t_axis','dtrace','dsg','fps','p2u', ...
               'halfwin','events'};
        % which = this quantity, with its axis label. The names are the settled
        % vocabulary; see BACKLOG.md
        quantities = {'divergence',      'divergence  \DeltaA/A'; ...
                      'strain_radial',   'strain_radial  n''En'; ...
                      'strain_hoop',     'strain_hoop  t''Et'; ...
                      'disp_radial',     'disp_radial  (+ outward, \mum)'; ...
                      'disp_tangential', 'disp_tangential  (+ CCW, \mum)'; ...
                      'volume_out',      'volume_out  (cumulative, \mum^2)'};
    end

    methods
        function obj = piv_figure(run, fields, profile, ctx)
            arguments
                run     struct
                fields  cell
                profile vfield_profile
                ctx     struct
            end
            if numel(run) ~= numel(fields)
                error('piv_figure:sizeMismatch', ...
                    'run has %d events but fields has %d entries.', numel(run), numel(fields));
            end
            if numel(run) ~= profile.n_row
                error('piv_figure:sizeMismatch', ...
                    'run has %d events but profile has %d rows.', numel(run), profile.n_row);
            end
            % row k must BE event k. append() fills in order, but a caller that
            % built the profile separately could have reordered it, and every
            % method below indexes both by the same k
            for k = 1:numel(run)
                if profile.rows(k).event_idx ~= run(k).id
                    error('piv_figure:rowMismatch', ...
                        'profile row %d is event %d but run(%d) is event %d.', ...
                        k, profile.rows(k).event_idx, k, run(k).id);
                end
            end
            missing = obj.req_fields(~isfield(ctx, obj.req_fields));
            if ~isempty(missing)
                error('piv_figure:ctxMissing', 'ctx is missing field(s): %s', ...
                    strjoin(missing, ', '));
            end
            obj.run = run;
            obj.fields = fields;
            obj.profile = profile;
            obj.ctx = ctx;

            obj.lims = struct();
            for q = 1:size(obj.quantities, 1)
                name = obj.quantities{q, 1};
                peak = nan(1, numel(run));
                for k = 1:numel(run)
                    v = profile.rows(k).cells.(name);
                    peak(k) = max(abs(v(:)), [], 'omitnan');
                end
                obj.lims.(name) = obj.symlim(peak);
            end
        end

        function q = quiveropts(obj)
            %QUIVEROPTS  Name-value cell for vfield_plotquiver. frame_span is pinned to 1.
            q = {"block_size",1, "frame_span",1, "scale",obj.scale, "linewidth",obj.linew, ...
                 "head_len",obj.headsz, "head_width",obj.headsz, ...
                 "filled_head",true, "head_abs",true};
        end

        function ax = vectors(obj, k, ax)
            %VECTORS  Ensemble displacement field over the FROM-state frame.
            if nargin < 3 || isempty(ax)
                figure('Name', obj.figname(k, 'piv'));
                ax = gca;
            end
            R = obj.run(k);
            imshow(obj.bg(k), 'Parent', ax);
            hold(ax, 'on');
            vfield_plotquiver(R.uv, obj.quiveropts{:});
            title(ax, obj.label(k));
        end

        function ax = map(obj, k, which, ax)
            %MAP  Overlay a per-cell quantity on the FROM frame, vectors on top.
            arguments
                obj
                k (1,1) double
                which (1,:) char = 'divergence'
                ax = []
            end
            if isempty(ax)
                figure('Name', obj.figname(k, which));
                ax = gca;
            end
            [M, cl, ttl] = obj.field(k, which);
            imshow(repmat(obj.bg(k), [1 1 3]), 'Parent', ax);
            hold(ax, 'on');
            h = imagesc(ax, M);
            set(h, 'AlphaData', obj.alpha * ~isnan(M));
            axis(ax, 'image');
            colormap(ax, obj.cmap);
            colorbar(ax);
            if isfinite(cl) && cl > 0
                clim(ax, [-cl cl]);
            end
            vfield_plotquiver(obj.run(k).uv, obj.quiveropts{:});
            title(ax, [ttl '   ' obj.label(k)]);
        end

        function ax = trace(obj, k, ax)
            %TRACE  Diameter around the event, with the merged window shaded.
            if nargin < 3 || isempty(ax)
                figure('Name', obj.figname(k, 'trace'));
                ax = gca;
            end
            R = obj.run(k);
            C = obj.ctx;
            w  = [C.events(R.id).start_f, C.events(R.id).end_f];
            rg = max(1, w(1) - 60):min(numel(C.dtrace), w(2) + 60);
            plot(ax, C.t_axis(rg), C.dtrace(rg), 'Color', [.75 .75 .75]);
            hold(ax, 'on');
            grid(ax, 'on');
            plot(ax, C.t_axis(rg), C.dsg(rg), 'k', 'LineWidth', 1.4);
            yl = ylim(ax);
            % Shade the merged bout window. plot_sleep_patches keys colour off the
            % field name, so this matches the review figure
            hp = plot_sleep_patches(ax, struct('trans_window', w), C.t_axis, ...
                'YLim', yl, 'FaceAlpha', .3);
            % Surrounding state bands, when ctx carries them
            if isfield(C, 'state_idx')
                st = {'roughnrem', 'rem', 'roughawake'};
                if isfield(C, 'states')
                    st = C.states;
                end
                bands = struct();
                for q = 1:numel(st)
                    if isfield(C.state_idx, st{q})
                        bands.(st{q}) = C.state_idx.(st{q});
                    end
                end
                if ~isempty(fieldnames(bands))
                    hp = [hp, plot_sleep_patches(ax, bands, C.t_axis, ...
                        'YLim', yl, 'FaceAlpha', .13)];
                end
            end
            uistack(hp, 'bottom');
            plot(ax, C.t_axis([R.from R.to]), C.dsg([R.from R.to]), 'g-', 'LineWidth', 2);
            plot(ax, C.t_axis(R.from), C.dsg(R.from), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 9);
            plot(ax, C.t_axis(R.to),   C.dsg(R.to),   'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 9);
            ylim(ax, yl);
            xlim(ax, C.t_axis([rg(1) rg(end)]));
            xlabel(ax, 't (s)');
            ylabel(ax, 'diameter (\mum)');
            title(ax, ['open=from, filled=to   ' obj.label(k)]);
        end

        function slices(obj, k)
            %SLICES  sliceViewer of the frames PIV ACTUALLY used, plus the slice map.
            %   Odd slices = FROM state, even = TO state, so stepping one slice at a
            %   time works as a blink comparator (wall motion vs whole-field drift).
            R = obj.run(k);
            C = obj.ctx;
            nT = numel(C.dtrace);
            wf = max(1, R.from - C.halfwin):min(nT, R.from + C.halfwin);
            wt = max(1, R.to   - C.halfwin):min(nT, R.to   + C.halfwin);
            kk = min(numel(wf), numel(wt));
            a  = round(linspace(wf(1), wf(end), kk));
            b  = round(linspace(wt(1), wt(end), kk));
            inter = zeros(size(C.S,1), size(C.S,2), 2*kk);
            inter(:,:,1:2:end) = C.S(:,:,a);
            inter(:,:,2:2:end) = C.S(:,:,b);
            fprintf('ev%-3d %-11s %-12s gap %5.1f s | FROM %s | TO %s\n', R.id, R.state, ...
                R.pol, R.rise_s, mat2str(a), mat2str(b));
            figure('Name', sprintf('ev%d %s %s  PIV frames (odd=from, even=to)', ...
                R.id, R.state, R.pol));
            sliceViewer(inter);
        end

        function panel(obj, k)
            %PANEL  disp_radial | divergence | volume_out | diameter, one figure.
            figure('Name', obj.figname(k, 'panel'), 'Position', [40 80 1750 430]);
            tl = tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
            obj.map(k, 'disp_radial', nexttile);
            obj.map(k, 'divergence',  nexttile);
            obj.map(k, 'volume_out',  nexttile);
            obj.trace(k,              nexttile);
            title(tl, obj.label(k), 'FontWeight', 'bold');
        end

        function tiles(obj, pol, which)
            %TILES  One column per event of a polarity; one row per quantity + trace.
            arguments
                obj
                pol   (1,:) char {mustBeMember(pol, {'dilation','constriction','none'})} = 'dilation'
                which cell = {'disp_radial','divergence','volume_out'}
            end
            sel = find(strcmp({obj.run.pol}, pol));
            if isempty(sel)
                warning('piv_figure:noEvents', 'no %s events', pol);
                return
            end
            nr = numel(which) + 1;
            figure('Name', sprintf('all %s events', pol), 'Position', [20 25 1860 980]);
            tl = tiledlayout(nr, numel(sel), 'TileSpacing', 'compact', 'Padding', 'compact');
            for row = 1:nr
                for ii = 1:numel(sel)
                    ax = nexttile;
                    if row <= numel(which)
                        obj.map(sel(ii), which{row}, ax);
                    else
                        obj.trace(sel(ii), ax);
                    end
                    if ii ~= numel(sel) || row > numel(which)
                        colorbar(ax, 'off');
                    end
                    if row > 1
                        title(ax, '');
                    end
                    set(ax, 'FontSize', 7);
                end
            end
            title(tl, sprintf('%s  (n=%d)   shared colour limits', upper(pol), numel(sel)), ...
                'FontWeight', 'bold');
        end

        function [bin_range, wedge_ok] = radprofile(obj, opt)
            %RADPROFILE  Per-event radial profile, on ONE wedge set for every row.
            %   The wedge set comes from vfield_profile, which requires a wedge to
            %   span the whole interval before it may contribute. Without that this
            %   plot averages more wedges near the wall than far from it and the
            %   slope it shows is partly the change in membership.
            arguments
                obj
                opt.quantity (1,:) char = 'volume_out'
                opt.bin_range (1,2) double = [NaN NaN]
            end
            allrows = true(1, obj.profile.n_row);
            bin_range = opt.bin_range;
            if any(isnan(bin_range))
                [bin_range, wedge_ok] = obj.profile.widest_common(allrows, opt.quantity);
            else
                wedge_ok = obj.profile.wedges_over(allrows, bin_range, opt.quantity);
            end
            if any(isnan(bin_range))
                warning('piv_figure:noCommonSupport', ...
                    'no wedge spans any interval for every row; nothing to plot');
                return
            end
            curves = obj.profile.curve_rows(allrows, bin_range, quantity = opt.quantity);
            r = obj.profile.bin_center_um;

            figure('Name', 'radial profile', 'Position', [120 80 1000 460]);
            hold on
            grid on
            for k = 1:numel(obj.run)
                switch obj.run(k).pol
                    case 'dilation';     co = [.85 .25 .25];
                    case 'constriction'; co = [.20 .40 .85];
                    otherwise;           co = [.55 .55 .55];
                end
                plot(r, curves(k, :), '-o', 'Color', co, 'LineWidth', 1.1, ...
                    'MarkerSize', 4, 'HandleVisibility', 'off');
            end
            yline(0, 'k--', 'HandleVisibility', 'off');
            xlabel('\mum from the vessel wall');
            ylabel(obj.axis_label(opt.quantity));
            title(sprintf(['red = dilation, blue = constriction, grey = quiet   |   ' ...
                '%d/%d wedges over %.1f-%.1f \\mum'], nnz(wedge_ok), obj.profile.n_wedge, ...
                r(bin_range(1)), r(bin_range(2))));
        end

        function signcheck(obj)
            %SIGNCHECK  dD vs measured direction; the two polarities must separate.
            %   err  ~isd used to mean constriction. A stable control (pol 'none')
            %        has no direction to be consistent WITH, so counting it as a
            %        failed constriction turned the controls into a fake defect
            %        rate. They are drawn, greyed, and left out of the count
            isd = strcmp({obj.run.pol}, 'dilation');
            isc = strcmp({obj.run.pol}, 'constriction');
            isn = ~isd & ~isc;
            mdr = obj.cell_mean('disp_radial');
            mdv = obj.cell_mean('divergence');
            figure('Name', 'sign check', 'Position', [200 200 950 400]);
            tiledlayout(1, 2, 'TileSpacing', 'compact');
            for p = 1:2
                if p == 1
                    % 1. Wall motion: dilation pushes out, constriction pulls in
                    y = mdr;
                    yl = 'mean disp\_radial (\mum, + outward)';
                    t = 'wall motion follows \DeltaD';
                    ok = (isd & y > 0) | (isc & y < 0);
                else
                    % 2. Volume: the sign the tissue's own divergence takes, which
                    %    is a RESULT and not a check -- see DIVERGENCE_RESULT.md.
                    %    Consistency here is counted, not asserted
                    y = mdv;
                    yl = 'mean divergence';
                    t = 'divergence follows polarity';
                    ok = (isd & y > 0) | (isc & y < 0);
                end
                nexttile;
                hold on
                grid on
                scatter([obj.run(isd).diameter_change], y(isd), 70, 'r', 'filled');
                scatter([obj.run(isc).diameter_change], y(isc), 70, 'b', 'filled');
                if any(isn)
                    scatter([obj.run(isn).diameter_change], y(isn), 70, ...
                        [.55 .55 .55], 'filled');
                end
                xline(0, 'k--');
                yline(0, 'k--');
                xlabel('\DeltaD (\mum)');
                ylabel(yl);
                % the denominator is the rows that HAVE a prescribed direction
                title(sprintf('%s   (%d/%d consistent)', t, nnz(ok), nnz(isd | isc)));
                if p == 1
                    leg = [{'dilation', 'constriction'}, ...
                           repmat({'stable control'}, 1, any(isn))];
                    legend(leg, 'Location', 'best');
                end
            end
        end

        function T = summary(obj)
            %SUMMARY  One row per event; returns the table as well as printing it.
            n = numel(obj.run);
            T = table('Size', [n 9], 'VariableTypes', ...
                {'double','string','string','double','double','double','double','double','double'}, ...
                'VariableNames', {'id','state','pol','dur_s','dD_um', ...
                                  'disp_radial','divergence','frac_neg','n_wedge_data'});
            mdr = obj.cell_mean('disp_radial');
            mdv = obj.cell_mean('divergence');
            for k = 1:n
                cells = obj.profile.rows(k).cells;
                v = cells.divergence(~isnan(cells.divergence));
                % wedges holding ANY data. NOT the number radprofile draws with,
                % which is the smaller set that spans the whole interval. The two
                % being confused is the bug this class was rebuilt around
                with_data = nnz(isfinite(cells.reach_bin));
                T(k,:) = {obj.run(k).id, string(obj.run(k).state), string(obj.run(k).pol), ...
                    obj.run(k).rise_s, obj.run(k).diameter_change, ...
                    mdr(k), mdv(k), mean(v < 0), with_data};
            end
            disp(T);
        end
    end

    methods (Access = private)
        function s = label(obj, k)
            %LABEL  Event caption. TeX stays outside the sprintf format string.
            R = obj.run(k);
            s = [sprintf('ev%d  %s bout %d  %s  %.1f s  ', R.id, R.state, R.bout, ...
                 R.pol, R.rise_s) '\DeltaD=' sprintf('%+.2f', R.diameter_change) ' \mum'];
        end

        function s = figname(obj, k, tag)
            R = obj.run(k);
            s = sprintf('%s ev%d %s %s', tag, R.id, R.state, R.pol);
        end

        function b = bg(obj, k)
            b = mat2gray(obj.ctx.S(:, :, obj.run(k).from));   % the EARLIER state
        end

        function [M, cl, ttl] = field(obj, k, which)
            % paint() spreads the per-cell value over the frame. Display only:
            % vfield_polar does not keep a dense map and nothing reads one back
            ttl = obj.axis_label(which);
            cells = obj.profile.rows(k).cells;
            M = obj.fields{k}.paint(cells.(which));
            cl = obj.lims.(which);
        end

        function m = cell_mean(obj, quantity)
            % One number per event: the mean over every cell that has a value.
            % Not a radial profile -- for a profile use radprofile(), which fixes
            % the wedge set. Here every cell of one event is one sample and the
            % membership question does not arise
            m = nan(1, numel(obj.run));
            for k = 1:numel(obj.run)
                v = obj.profile.rows(k).cells.(quantity);
                m(k) = mean(v(:), 'omitnan');
            end
        end

        function s = axis_label(obj, which)
            hit = strcmp(obj.quantities(:,1), which);
            if ~any(hit)
                error('piv_figure:unknownQuantity', ...
                    '%s is not one of: %s', which, strjoin(obj.quantities(:,1)', ', '));
            end
            s = obj.quantities{hit, 2};
        end
    end

    methods (Static, Access = private)
        function v = symlim(x)
            v = max(x(isfinite(x)));
            if isempty(v) || v <= 0
                v = 1;
            end
        end
    end
end
