function plot_radon_centershift(t_axis,centershift_amp,pixel2um,clee)
    centershift_fig = make_fig("Center shift","normal");
    centershift_fig.resolution = pixel2um;
    centershift_fig.update_figsize([9 3]);
    centershift_fig.loc.x = t_axis;
    centershift_fig.plot_line(centershift_amp,clee.lch(45,80,5));
end