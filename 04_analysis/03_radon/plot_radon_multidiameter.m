function radon_multidia = plot_radon_multidiameter(t_axis,diameters,var_diameter,pixel2um,clee)       
    radon_multidia = make_fig("multiple diameter");
    radon_multidia.update_figsize([9 3]); % 6:1 ratio for 2 subplots (each 3:1)
    radon_multidia.loc.x = t_axis;
    radon_multidia.resolution = pixel2um;
    hold(radon_multidia.ax,"on")
    [~,maxvar_angle] = max(var_diameter);
    %%

    %%
    radon_multidia.reset_axis
    tmp.init_val = 0;
    for angleidx = 0:5
        tmp.plot_angle = mod(maxvar_angle+angleidx*30,180);
        tmp.hue_angle = angleidx*50+5;
        tmp.init_val = 30*angleidx - diameters(tmp.plot_angle,1); 
        tmp.plot_data = diameters(tmp.plot_angle,:) +tmp.init_val;
        plot(radon_multidia.ax,t_axis, tmp.plot_data.*pixel2um , color = clee.lch(45,80,tmp.hue_angle))
        text(radon_multidia.ax,-10,(tmp.plot_data(1)+5).*pixel2um ,strcat('Angle',string(tmp.plot_angle)), "HorizontalAlignment",'right')
    end
    %%
    radon_multidia.change_xylim([-80 max(t_axis)],...
        [(min(diameters(maxvar_angle,:))- diameters(maxvar_angle,1)).*pixel2um,max(tmp.plot_data.*pixel2um)])
    
    set(radon_multidia.ax,"YTick",[])
    radon_multidia.ax.YAxis.Visible = 'off';
    radon_multidia.put_xaxistitle("Time (sec)")
end

