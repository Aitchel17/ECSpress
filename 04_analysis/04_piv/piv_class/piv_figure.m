classdef piv_figure < handle
%PIV_FIGURE  Display layer for per-event ensemble-PIV results.
%   Holds the two settings that must stay constant across events:
%     quiver scale   arrow length = displacement(px) * frame_span * scale, and
%                    frame_span is pinned to 1 because uv is already a total
%                    from->to displacement (see piv_run_events)
%     colour limits  computed once over all events, symmetric about 0, so the maps
%                    show amplitude differences instead of self-normalising
%
% IN   run    1xN struct from piv_run_events
%             (id,state,pol,from,to,rise_s,diameter_change,nfr,uv,corr)
%      polar  1xN cell from piv_polar_events
%             (cells,rcen,vr,vt,pdiv,vrmap,vtmap,div,info)
%      ctx    struct: S, coremask, exclmask, t_axis, dtrace, dsg, fps, p2u,
%             halfwin, events, and optionally state_idx / states
% OUT  obj    handle object; call the methods below
%
%   USAGE
%     F = piv_figure(piv_run, polar, ctx);
%     F.summary();          % text table
%     F.map(2, 'radial');   % dv_r/dr overlay + quiver for event index 2
%     F.panel(2);           % vr | dv_r/dr | plain div | diameter, one figure
%     F.tiles('dilation');  % every dilation event tiled
%     F.radprofile();  F.signcheck();  F.slices(2);

    properties
        run                     % per-event PIV results
        polar                   % per-event polar/divergence maps
        ctx                     % shared context (stack, masks, traces, geometry)
        lims                    % struct(vr,dvr,div): symmetric limits shared by all events
        scale   = 5             % quiver length multiplier (frame_span is always 1)
        alpha   = 0.55          % overlay transparency
        cmap    = 'turbo'
        headsz  = 5             % absolute arrowhead size (px)
        linew   = 1
    end

    properties (Constant, Access = private)
        req_fields = {'S','coremask','exclmask','t_axis','dtrace','dsg','fps','p2u', ...
               'halfwin','events'};
    end

    methods
        function obj = piv_figure(run, polar, ctx)
            arguments
                run   struct
                polar cell
                ctx   struct
            end
            if numel(run) ~= numel(polar)
                error('piv_figure:sizeMismatch', ...
                    'run has %d events but polar has %d entries.', numel(run), numel(polar));
            end
            missing = obj.req_fields(~isfield(ctx, obj.req_fields));
            if ~isempty(missing)
                error('piv_figure:ctxMissing', 'ctx is missing field(s): %s', ...
                    strjoin(missing, ', '));
            end
            obj.run   = run;
            obj.polar = polar;
            obj.ctx   = ctx;
            obj.lims  = struct( ...
                'vr',  obj.symlim(cellfun(@(A) max(abs(A.vrmap(:)), [], 'omitnan'), polar)), ...
                'dvr', obj.symlim(cellfun(@(A) max(abs(A.pdiv(:)),  [], 'omitnan'), polar)), ...
                'div', obj.symlim(cellfun(@(A) prctile(abs(A.div(~isnan(A.div))), 99), polar)), ...
                'vt',  obj.symlim(cellfun(@(A) max(abs(A.vtmap(:)), [], 'omitnan'), polar)));
        end

        function q = quiveropts(obj)
            %QUIVEROPTS  Name-value cell for vfield_plotquiver. frame_span is pinned to 1.
            q = {"block_size",1, "frame_span",1, "scale",obj.scale, "linewidth",obj.linew, ...
                 "head_len",obj.headsz, "head_width",obj.headsz, ...
                 "filled_head",true, "head_abs",true};
        end

        function ax = vectors(obj, k, ax)
            %VECTORS  Ensemble displacement field over the FROM-state frame.
            if nargin < 3 || isempty(ax); figure('Name', obj.figname(k, 'piv')); ax = gca; end
            R = obj.run(k);
            imshow(obj.bg(k), 'Parent', ax); hold(ax, 'on');
            vfield_plotquiver(R.uv, obj.quiveropts{:});
            title(ax, obj.label(k));
        end

        function ax = map(obj, k, which, ax)
            %MAP  Overlay a scalar field on the FROM frame and draw the vectors on top.
            %   which: 'radial' = dv_r/dr (- compressed), 'plain' = du/dx+dv/dy,
            %          'vr' = mean radial velocity (+ outward),
            %          'vt' = mean tangential velocity (+ counter-clockwise).
            arguments
                obj, k (1,1) double
                which (1,:) char {mustBeMember(which, {'radial','plain','vr','vt'})} = 'radial'
                ax = []
            end
            if isempty(ax); figure('Name', obj.figname(k, which)); ax = gca; end
            [M, cl, ttl] = obj.field(k, which);
            imshow(repmat(obj.bg(k), [1 1 3]), 'Parent', ax); hold(ax, 'on');
            h = imagesc(ax, M); set(h, 'AlphaData', obj.alpha * ~isnan(M));
            axis(ax, 'image'); colormap(ax, obj.cmap); colorbar(ax);
            if isfinite(cl) && cl > 0; clim(ax, [-cl cl]); end
            vfield_plotquiver(obj.run(k).uv, obj.quiveropts{:});
            title(ax, [ttl '   ' obj.label(k)]);
        end

        function ax = trace(obj, k, ax)
            %TRACE  Diameter around the event, with the merged window shaded.
            if nargin < 3 || isempty(ax); figure('Name', obj.figname(k, 'trace')); ax = gca; end
            R = obj.run(k);  C = obj.ctx;
            w  = [C.events(R.id).start_f, C.events(R.id).end_f];
            rg = max(1, w(1) - 60):min(numel(C.dtrace), w(2) + 60);
            plot(ax, C.t_axis(rg), C.dtrace(rg), 'Color', [.75 .75 .75]); hold(ax, 'on'); grid(ax, 'on');
            plot(ax, C.t_axis(rg), C.dsg(rg), 'k', 'LineWidth', 1.4);
            yl = ylim(ax);
            % Shade the merged bout window. plot_sleep_patches keys colour off the
            % field name, so this matches the review figure
            hp = plot_sleep_patches(ax, struct('trans_window', w), C.t_axis, ...
                'YLim', yl, 'FaceAlpha', .3);
            % Surrounding state bands, when ctx carries them
            if isfield(C, 'state_idx')
                st = {'roughnrem', 'rem', 'roughawake'};
                if isfield(C, 'states'); st = C.states; end
                bands = struct();
                for q = 1:numel(st)
                    if isfield(C.state_idx, st{q}); bands.(st{q}) = C.state_idx.(st{q}); end
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
            ylim(ax, yl); xlim(ax, C.t_axis([rg(1) rg(end)]));
            xlabel(ax, 't (s)'); ylabel(ax, 'diameter (\mum)');
            title(ax, ['open=from, filled=to   ' obj.label(k)]);
        end

        function slices(obj, k)
            %SLICES  sliceViewer of the frames PIV ACTUALLY used, plus the slice map.
            %   Odd slices = FROM state, even = TO state, so stepping one slice at a
            %   time works as a blink comparator (wall motion vs whole-field drift).
            R = obj.run(k);  C = obj.ctx;  nT = numel(C.dtrace);
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
            %PANEL  vr | dv_r/dr | plain divergence | diameter, one figure.
            figure('Name', obj.figname(k, 'panel'), 'Position', [40 80 1750 430]);
            tl = tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
            obj.map(k, 'vr',     nexttile);
            obj.map(k, 'radial', nexttile);
            obj.map(k, 'plain',  nexttile);
            obj.trace(k,         nexttile);
            title(tl, obj.label(k), 'FontWeight', 'bold');
        end

        function tiles(obj, pol, which)
            %TILES  One column per event of a polarity; rows = vr / dv_r/dr / div / trace.
            arguments
                obj
                pol   (1,:) char {mustBeMember(pol, {'dilation','constriction','none'})} = 'dilation'
                which cell = {'vr','radial','plain'}
            end
            sel = find(strcmp({obj.run.pol}, pol));
            if isempty(sel); warning('piv_figure:noEvents', 'no %s events', pol); return; end
            nr = numel(which) + 1;
            figure('Name', sprintf('all %s events', pol), 'Position', [20 25 1860 980]);
            tl = tiledlayout(nr, numel(sel), 'TileSpacing', 'compact', 'Padding', 'compact');
            for row = 1:nr
                for ii = 1:numel(sel)
                    ax = nexttile;
                    if row <= numel(which); obj.map(sel(ii), which{row}, ax);
                    else;                   obj.trace(sel(ii), ax);
                    end
                    if ii ~= numel(sel) || row > numel(which); colorbar(ax, 'off'); end
                    if row > 1; title(ax, ''); end
                    set(ax, 'FontSize', 7);
                end
            end
            title(tl, sprintf('%s  (n=%d)   shared colour limits', upper(pol), numel(sel)), ...
                'FontWeight', 'bold');
        end

        function radprofile(obj)
            %RADPROFILE  dv_r/dr vs distance from the vessel wall, all events.
            figure('Name', 'radial profile', 'Position', [120 80 1000 460]); hold on; grid on
            for k = 1:numel(obj.run)
                A = obj.polar{k};
                % err  a two-way test coloured the stable controls as constrictions
                switch obj.run(k).pol
                    case 'dilation';     co = [.85 .25 .25];
                    case 'constriction'; co = [.20 .40 .85];
                    otherwise;           co = [.55 .55 .55];
                end
                plot(mean(A.rcen, 2, 'omitnan'), mean(A.cells, 2, 'omitnan'), '-o', ...
                    'Color', co, 'LineWidth', 1.1, 'MarkerSize', 4, 'HandleVisibility', 'off');
            end
            yline(0, 'k--', 'HandleVisibility', 'off');
            xlabel('distance from vessel wall (px)'); ylabel('dv_r/dr');
            title('red = dilation (PVS compressed), blue = constriction (refill), grey = stable control');
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
            mvr = cellfun(@(A) mean(A.vr(:),    'omitnan'), obj.polar);
            mdv = cellfun(@(A) mean(A.cells(:), 'omitnan'), obj.polar);
            figure('Name', 'sign check', 'Position', [200 200 950 400]);
            tiledlayout(1, 2, 'TileSpacing', 'compact');
            for p = 1:2
                if p == 1
                    % 1. Wall motion: dilation pushes out, constriction pulls in
                    y = mvr;  yl = 'mean v_r (px, + outward)';  t = 'wall motion follows \DeltaD';
                    ok = (isd & y > 0) | (isc & y < 0);
                else
                    % 2. PVS strain: v_r decays with distance, so dilation squeezes
                    y = mdv;  yl = 'mean dv_r/dr';  t = 'PVS compressed on dilation';
                    ok = (isd & y < 0) | (isc & y > 0);
                end
                nexttile; hold on; grid on
                scatter([obj.run(isd).diameter_change], y(isd), 70, 'r', 'filled');
                scatter([obj.run(isc).diameter_change], y(isc), 70, 'b', 'filled');
                if any(isn)
                    scatter([obj.run(isn).diameter_change], y(isn), 70, ...
                        [.55 .55 .55], 'filled');
                end
                xline(0, 'k--'); yline(0, 'k--');
                xlabel('\DeltaD (\mum)'); ylabel(yl);
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
            T = table('Size', [n 8], 'VariableTypes', ...
                {'double','string','string','double','double','double','double','double'}, ...
                'VariableNames', {'id','state','pol','dur_s','dD_um','v_r_px','dvr_dr','frac_neg'});
            for k = 1:n
                A = obj.polar{k};  v = A.cells(~isnan(A.cells));
                T(k,:) = {obj.run(k).id, string(obj.run(k).state), string(obj.run(k).pol), ...
                    obj.run(k).rise_s, obj.run(k).diameter_change, mean(A.vr(:), 'omitnan'), ...
                    mean(v), mean(v < 0)};
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
            A = obj.polar{k};
            switch which
                case 'radial'; M = A.pdiv;  cl = obj.lims.dvr; ttl = 'dv_r/dr  (- compressed)';
                case 'plain';  M = A.div;   cl = obj.lims.div; ttl = 'div  \DeltaA/A';
                case 'vr';     M = A.vrmap; cl = obj.lims.vr;  ttl = 'v_r  (+ outward, px)';
                case 'vt';     M = A.vtmap; cl = obj.lims.vr;  ttl = ['v_' char(92) 'theta  (+ CCW, px)'];
            end
        end
    end

    methods (Static, Access = private)
        function v = symlim(x)
            v = max(x(isfinite(x)));
            if isempty(v) || v <= 0; v = 1; end
        end
    end
end
