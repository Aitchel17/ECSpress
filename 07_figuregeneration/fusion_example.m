%FUSION_EXAMPLE  Building one composite by hand, in the order the methods run.
%   Nothing here is wrapped. Each line is one method call and the order they appear
%   in IS the order they have to happen in: the board exists before a panel can go
%   on it, every panel is on before a column can be linked, and the link is made
%   before the file is written or the bottom panel keeps its ticks.
%
%   Run review_session first. That is what leaves pax_fig, figStruct and sleepscore
%   in the workspace, and this reads them straight out of it.

%% 1. the panels, each in its own window, exactly as the makefigs leave them
sleepstrip = plot_sleep_strip(sleepscore);
pax_fig = analysis_pax_makefig(session.pax_fwhm, session.pax_fwhm.t_axis, ...
    session.img_param.pixel2um, session.dir_struct.figures_fwhm);

%% 2. what is open and can be snapped on
fusion_figure.available

%% 3. a board, then one panel per cell
% grid is [rows columns]. Left column carries time, right column carries the panels
% that have no abscissa, which is why link_x is called on column 1 and not on 2.
board = fusion_figure('grid', [3 2], 'width', 10, 'height', 9, 'fontsize', 9);
board.add(sleepstrip,               'at', [1 1], 'label', "A");
board.add(figStruct.spectrogram.ax, 'at', [2 1], 'label', "B", 'colorbar', true);
board.add(pax_fig.pvs.ax,           'at', [3 1], 'label', "C");
board.link_x(1);

%% 4. the right column, one panel at a time
board.add(pax_fig.bv.ax, 'at', [1 2], 'label', "D");

%% 5. wrong panel, take it off and put another in its place
board.drop([1 2]);
board.add(pax_fig.totalpvs_changes.ax, 'at', [1 2], 'label', "D");

%% 6. write it
dirs = dirs_central();
board.save_figure(fullfile(dirs.secondary_root, "merged_igkl_igkltdt"), "fusion_example");
