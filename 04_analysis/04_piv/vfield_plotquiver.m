function vfield_plotquiver(flow, opt)
%VFIELD_PLOTQUIVER  PIVlab-style quiver plot of a displacement field.
%   4D input is summed over the pair axis (time_range selects which pairs); a 2D
%   ensemble field carries per-frame displacement, so frame_span scales it to the
%   total over a window. Leave frame_span at 1 for 4D input. Draws into gca.
%
% IN   flow         HxWxNx2 pair stack, or HxWx1x2 / HxWx2 single field
%      time_range   [start end] pair indices to sum, 4D only, default all pairs
%      frame_span   frames the field is multiplied by, 2D only, default 1
%      scale        vector scale factor, default 1
%      block_size   spatial averaging block, default 1 = no averaging
%      uniform      normalize to unit length, default false
%      color_valid  RGB quiver colour, default [0 1 0]
%      linewidth    arrow / shaft line width, default 0.5
%      max_head     arrowhead size rel. to arrow length, quiver path only,
%                   default 0.2, 0 = plain line segments
%      filled_head  thin shaft plus filled triangular head, default false
%      head_abs     head sizes are absolute px, not a fraction of arrow length,
%                   default false
%      head_len     filled-head length, fraction of arrow length or px, default 0.35
%      head_width   filled-head base width, fraction of arrow length or px, default 0.22

arguments
    flow             {mustBeNumeric}
    opt.time_range   double  = []   % [start end] pairs to sum (4D only)
    opt.frame_span   double  = 1    % 2D ensemble: multiply per-frame field by #frames
    opt.scale        double  = 1
    opt.block_size   double  = 1
    opt.uniform      logical = false
    opt.color_valid  double  = [0 1 0]
    opt.linewidth    double  = 0.5
    opt.max_head     double  = 0.2   % quiver-path arrowhead size rel. to length (0 = no heads)
    opt.filled_head  logical = false % true = thin shaft + filled triangle head (decoupled)
    opt.head_len     double  = 0.35  % filled-head length  (fraction of arrow len, or px if head_abs)
    opt.head_width   double  = 0.22  % filled-head base width (fraction of arrow len, or px if head_abs)
    opt.head_abs     logical = false % true = uniform absolute head size (px) for every arrow
end

% 1. Reduce the input to one HxW displacement field
ndims_flow = ndims(flow);
if ndims_flow == 4
    % 1.1 Pair indices to accumulate, clamped to the stack
    N = size(flow, 3);
    if isempty(opt.time_range)
        t_idx = 1:N;
    else
        t_idx = round(opt.time_range(1)) : round(opt.time_range(2));
        t_idx = t_idx(t_idx >= 1 & t_idx <= N);
    end
    % 1.2 Sum displacement over those pairs
    u = sum(flow(:,:,t_idx,1), 3);
    v = sum(flow(:,:,t_idx,2), 3);
else
    % 1.3 Single or ensemble field: split the two components
    f = squeeze(flow);
    u = f(:,:,1);
    v = f(:,:,2);
    % 1.4 On a 2D field time_range only sets the window length (frames)
    if ~isempty(opt.time_range)
        opt.frame_span = opt.time_range(2) - opt.time_range(1);
    end
end

% 2. Per-frame displacement to the total over the frame window
u = u * opt.frame_span;
v = v * opt.frame_span;

% 3. Pixel grid of the field
[H, W] = size(u);
[x, y] = meshgrid(1:W, 1:H);

% 4. Normalize to unit length, keeping direction only
if opt.uniform
    mag = sqrt(u.^2 + v.^2);
    mag(mag == 0) = 1;
    u = u ./ mag;
    v = v ./ mag;
end

% 5. Spatial block averaging
%   NaN-safe mean over block_size x block_size tiles; the arrow sits at the tile centre.
s = round(opt.block_size);
if s > 1
    % 5.1 Trim to an integer multiple of the block size
    Hr = floor(H/s)*s;
    Wr = floor(W/s)*s;
    u  = u(1:Hr, 1:Wr);
    v  = v(1:Hr, 1:Wr);

    % 5.2 Average down rows, then across columns
    u = squeeze(mean(reshape(u, s, Hr/s, Wr), 1, 'omitnan'));
    v = squeeze(mean(reshape(v, s, Hr/s, Wr), 1, 'omitnan'));
    u = squeeze(mean(reshape(u', s, Wr/s, Hr/s), 1, 'omitnan'))';
    v = squeeze(mean(reshape(v', s, Wr/s, Hr/s), 1, 'omitnan'))';

    % 5.3 Move the grid to the centre of each block
    xc = (0.5:Wr/s) * s;
    yc = (0.5:Hr/s) * s;
    [x, y] = meshgrid(xc, yc);
end

% 6. Draw the arrows
hold on;
if opt.filled_head
    % 6.1 Shaft thickness decoupled from head size
    draw_filled_arrows(x, y, u*opt.scale, v*opt.scale, ...
        opt.color_valid, opt.linewidth, opt.head_len, opt.head_width, opt.head_abs);
else
    % 6.2 Built-in quiver, one LineWidth for shaft and head
    quiver(x, y, u*opt.scale, v*opt.scale, ...
        'Color', opt.color_valid, ...
        'AutoScale', 'off', ...
        'MaxHeadSize', opt.max_head, ...
        'LineWidth', opt.linewidth);
end
hold off;
axis image;

end

% ---------------------------------------------------------------------------
function draw_filled_arrows(x, y, u, v, col, lw, head_len, head_width, head_abs)
%DRAW_FILLED_ARROWS  Thin-shaft arrows with filled triangular heads.
%   head_len/head_width are fractions of each arrow's length, or absolute pixels
%   when head_abs. NaN and zero-length vectors are skipped.

    % 1. Flatten and keep the drawable vectors
    x = x(:); y = y(:); u = u(:); v = v(:);
    L = hypot(u, v);
    ok = isfinite(L) & L > 0 & isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok); u = u(ok); v = v(ok); L = L(ok);
    if isempty(x); return; end

    % 2. Arrow frame and head geometry
    dx = u ./ L;  dy = v ./ L;              % unit along-arrow
    px = -dy;     py = dx;                  % unit perpendicular
    tipx = x + u; tipy = y + v;             % arrow tip
    if head_abs
        hl = min(head_len, L);              % uniform px length, clamped <= arrow length
        hw = head_width * ones(size(L));    % uniform px width
    else
        hl = head_len   .* L;               % head length     scales with arrow
        hw = head_width .* L;               % head base width  scales with arrow
    end
    bx = tipx - dx .* hl;  by = tipy - dy .* hl;   % head base centre

    % 3. Shafts: tail -> head base, NaN-separated into one line object
    n  = numel(x);
    sx = [x, bx, nan(n, 1)].';
    sy = [y, by, nan(n, 1)].';
    line(sx(:), sy(:), 'Color', col, 'LineWidth', lw);

    % 4. Heads: filled triangles [tip, left base, right base]
    lx = bx + px .* (hw / 2);  ly = by + py .* (hw / 2);
    rx = bx - px .* (hw / 2);  ry = by - py .* (hw / 2);
    V = [tipx tipy; lx ly; rx ry];                    % 3n x 2 vertices
    F = [(1:n).', (n + 1:2 * n).', (2 * n + 1:3 * n).'];  % n x 3 faces
    patch('Faces', F, 'Vertices', V, 'FaceColor', col, 'EdgeColor', 'none');
end
