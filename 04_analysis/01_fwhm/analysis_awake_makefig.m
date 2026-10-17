function fig = analysis_awake_makefig(time_table, t_axis, bv_px, totalpvs_px, pixel2um, ball_t, ball_velocity, save_dir, opt)
%ANALYSIS_AWAKE_MAKEFIG  The puff x locomotion bouts as a strip, and the ball velocity, vessel diameter and PVS thickness with the bouts behind them.
%   Caller: b_reviewsession
%
% IN   time_table     1x1 struct   bout tables in seconds, from bout_puffmove
%      t_axis         1 x T float  frame times of the FWHM series (s)
%      bv_px          1 x T float  vessel diameter, lower boundary minus upper (px)
%      totalpvs_px    1 x T float  PVS thickness, both sides (px)
%      pixel2um       1x1 float    um per pixel
%      ball_t         1 x B float  times of the ball series (s), analysis_analog ball.rs_taxis
%      ball_velocity  1 x B float  signed ball velocity (cm/s), analysis_analog ball.resampled_ball
%      save_dir       1x1 str      folder each window's svg is written to
%      opt.window_inch  1x2 float  size of the three trace windows (in)
%      opt.strip_inch   1x2 float  size of the strip window (in)
%      opt.bar_frac     1x1 float  puff bar height, as a fraction of the trace window's y range
% OUT  fig.strip    make_fig   a puff row and a locomotion / rest row, a legend of the seven categories
%      fig.ball     make_fig   ball velocity (cm/s)
%      fig.bv       make_fig   vessel diameter (um)
%      fig.pvs      make_fig   PVS thickness (um)
%   One window each, one axes each, one svg each: awakeStrip_figure,
%   ballVelocity_awake_figure, paxBV_awake_figure, paxPVS_awake_figure.
    arguments
        time_table (1,1) struct
        t_axis (1,:) double
        bv_px (1,:) double
        totalpvs_px (1,:) double
        pixel2um (1,1) double {mustBePositive}
        ball_t (1,:) double
        ball_velocity (1,:) double
        save_dir {mustBeTextScalar}
        opt.window_inch (1,2) double = [8 3]
        opt.strip_inch (1,2) double = [8 1.6]
        opt.bar_frac (1,1) double = 0.06
    end

    clee = color_lee;
    recording_end = t_axis(end);
    headroom_frac = 1.5 * opt.bar_frac;

    %% 1. The seven categories as two rows; locomotion and rest never overlap, so they share one
    fig.strip = make_fig('awakeStrip_figure');
    fig.strip.update_figsize(opt.strip_inch);
    fig.strip.reset_axis();
    patch_list = plot_awake_patches(fig.strip.ax, time_table, rows = true);
    xlim(fig.strip.ax, [0 recording_end]);
    ylim(fig.strip.ax, [0.4 2.6]);
    set(fig.strip.ax, 'YTick', 1:2, 'YTickLabel', {'locomotion / rest', 'puff'});
    fig.strip.put_xaxistitle('Time (sec)');
    key = legend(fig.strip.ax, patch_list, 'Orientation', 'horizontal', ...
        'Location', 'northoutside', 'Box', 'off', 'Interpreter', 'none');
    key.FontName = fig.strip.fontname;
    key.FontSize = fig.strip.fontsize;
    fig.strip.save2svg(save_dir);

    %% 2. Ball velocity
    ball_limit = padded_limits(ball_velocity, headroom_frac);
    fig.ball = make_fig('ballVelocity_awake_figure');
    fig.ball.update_figsize(opt.window_inch);
    fig.ball.reset_axis();
    fig.ball.loc.x = ball_t;
    fig.ball.plot_line(ball_velocity, clee.clist.black);
    xlim(fig.ball.ax, [0 recording_end]);
    ylim(fig.ball.ax, ball_limit);
    plot_awake_patches(fig.ball.ax, time_table, ylim = ball_limit, bar_frac = opt.bar_frac);
    fig.ball.put_yaxistitle('Ball velocity (cm/s)');
    fig.ball.put_xaxistitle('Time (sec)');
    fig.ball.save2svg(save_dir);

    %% 3. Vessel diameter
    bv_um = bv_px * pixel2um;
    bv_limit = padded_limits(bv_um, headroom_frac);
    fig.bv = make_fig('paxBV_awake_figure');
    fig.bv.update_figsize(opt.window_inch);
    fig.bv.resolution = pixel2um;
    fig.bv.reset_axis();
    fig.bv.loc.x = t_axis;
    fig.bv.plot_line(bv_px, clee.clist.magenta);
    ylim(fig.bv.ax, bv_limit);
    plot_awake_patches(fig.bv.ax, time_table, ylim = bv_limit, bar_frac = opt.bar_frac);
    fig.bv.put_yaxistitle('Vessel diameter (\mum)');
    fig.bv.put_xaxistitle('Time (sec)');
    fig.bv.save2svg(save_dir);

    %% 4. PVS thickness
    pvs_um = totalpvs_px * pixel2um;
    pvs_limit = padded_limits(pvs_um, headroom_frac);
    fig.pvs = make_fig('paxPVS_awake_figure');
    fig.pvs.update_figsize(opt.window_inch);
    fig.pvs.resolution = pixel2um;
    fig.pvs.reset_axis();
    fig.pvs.loc.x = t_axis;
    fig.pvs.plot_line(totalpvs_px, clee.clist.darkgreen);
    ylim(fig.pvs.ax, pvs_limit);
    plot_awake_patches(fig.pvs.ax, time_table, ylim = pvs_limit, bar_frac = opt.bar_frac);
    fig.pvs.put_yaxistitle('PVS thickness (\mum)');
    fig.pvs.put_xaxistitle('Time (sec)');
    fig.pvs.save2svg(save_dir);

    fprintf('Saved awakeStrip_figure, ballVelocity_awake_figure, paxBV_awake_figure, paxPVS_awake_figure to %s\n', save_dir);
end

function y_limit = padded_limits(value, headroom_frac)
%PADDED_LIMITS  [min max] of a series, with room above for the puff bars.
%   Caller: analysis_awake_makefig
%
% IN   value          1 x T float  the series
%      headroom_frac  1x1 float    extra height above the maximum, as a fraction of the range
% OUT  y_limit        1x2 float
    low = min(value);
    high = max(value);
    y_limit = [low, high + headroom_frac * (high - low)];
end
