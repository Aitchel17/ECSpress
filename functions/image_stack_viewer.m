function [roiCell, ax] = image_stack_viewer(img, varargin)
% 'colormap'           - Colormap name (char). Default: 'gray'.
%                         Supports special maps like 'icg', 'sO2'.
% 'scale'              - Intensity scaling: 'linear' (default), 'log', 
%                         'db', 'icg'.
% 'select'             - ROI selection mode: 'none' (default), 
%                         'rectangle', 'polygon', 'mask'.
% 'propertyName'       - Name of property being visualized (affects defaults).
% 'tag'                - Tag string for the ROI selection.
% 'plot'               - Logical (default: false). If true, plots time 
%                         series for selected ROI.
% 'filtering'          - Cell array {type, params} for image filtering.
% 'colorLimits'        - [min max] limits for color scaling.
% 'mask'               - Binary mask or numeric mask to overlay.
% 'propertyTitle'      - Title of the figure window.
% 'axesCoords'         - Cell array {x_axis_vec, z_axis_vec}.
% 'xlabel', 'zlabel'   - Labels for the axes.
% 'axesTitle'          - Title for the axes.
% 'frameTitles'        - Cell array of strings, one per frame.
% 'depthScatterNorm'   - Logical, apply depth scatter normalization.
% 'backNorm'           - Logical, apply background normalization.
% 'stimulationIndices' - Logical array indicating stimulation frames (for title color).
% 'stimulationTimeIndices' - Time values for individual block for frame analysis
% 'svdFilter'          - SVD filtering range for spatiotemporal filtering
% 'hardColorLimits'    - [min max] limits for color scaling that are not changed by user interaction
% 'printStats'         - Print statistics for the 3D tensor
% 'locCount'           - Number of Points to Localize in Each Frame
% 'fwhm'               - FWHM of Bubbles for Localization
% 'nLocalMax'          - Number of Local Maxima for Localization

  p = inputParser;
  addParameter(p, 'colormap', 'gray', @ischar);
  addParameter(p, 'scale', 'linear', @(x) ismember(x, {'linear', 'log', 'db', 'icg', 'exp'}));
  addParameter(p, 'expScale', 1, @(x) isnumeric(x));
  addParameter(p, 'select', 'none', @(x) ismember(x, {'none', 'rectangle', 'polygon', 'mask'}));
  addParameter(p, 'propertyName', '', @ischar);
  addParameter(p, 'tag', 'none', @ischar);
  addParameter(p, 'plot', false, @islogical);
  addParameter(p, 'filtering', {0, 0}, @(x) iscell(x) || isnumeric(x));
  addParameter(p, 'colorLimits', [], @(x) isnumeric(x) && numel(x) == 2);
  addParameter(p, 'mask', [], @(x) isnumeric(x) || islogical(x));
  addParameter(p, 'propertyTitle', 'Image Stack Viewer', @ischar);
  addParameter(p, 'axesCoords', {}, @iscell);
  addParameter(p, 'xlabel', '', @ischar);
  addParameter(p, 'zlabel', '', @ischar);
  addParameter(p, 'axesTitle', '', @ischar);
  addParameter(p, 'frameTitles', {}, @iscell);
  addParameter(p, 'depthScatterNorm', false, @islogical);
  addParameter(p, 'backNorm', false, @islogical);
  addParameter(p, 'stimulationIndices', [], @(x) isnumeric(x) || islogical(x));
  addParameter(p, 'stimulationTimeIndices', [], @(x) isnumeric(x));
  addParameter(p, 'svdFilter', [], @(x) isnumeric(x));
  addParameter(p, 'hardColorLimits', [], @(x) isnumeric(x) && numel(x) == 2);
  addParameter(p, 'printStats', true, @islogical);
  addParameter(p, 'locCount', 0, @(x) isnumeric(x));
  addParameter(p, 'fwhm', 3, @(x) isnumeric(x));
  addParameter(p, 'nLocalMax', 3, @(x) isnumeric(x));
  addParameter(p, 'locParams', {'g', 8}, @(x) iscell(x) && numel(x) == 2);
  parse(p, varargin{:});
  axesCoords = p.Results.axesCoords;
  cmap = p.Results.colormap;
  if startsWith(cmap, 'slanCM')
    if length(cmap) > 7
      cmapName = cmap(8:end);
    else
      cmapName = 'viridis';
    end
    
    flipMap = false;
    if length(cmapName) > 2 && strcmp(cmapName(end-1:end), '_r')
      cmapName = cmapName(1:end-2);
      flipMap = true;
    end
    
    try
      cmap = slanCM(cmapName, 256);
      if flipMap
        cmap = flipud(cmap);
      end
    catch ME
      warning('Failed to Load slanCM Colormap: %s. Using Default Gray.', ME.message);
      cmap = gray(256);
    end
  else
    if isequal(cmap, 'icg')
      cmap = icgColorMap();
    elseif isequal(cmap, 'sO2')
      cmap = sO2ColorMap();
    end
  end
  scale = p.Results.scale;
  select = p.Results.select;
  expScale = p.Results.expScale;
  roi = [];
  tag = p.Results.tag;
  seriesPlot = p.Results.plot;
  if seriesPlot && strcmp(select, 'none')
    warning('No ROI Selection Method Chosen. Proceeding Without ROI Selection.');
  end
  propertyName = p.Results.propertyName;
  propertyTitle = p.Results.propertyTitle;
  frameTitles = p.Results.frameTitles;
  stimulationIndices = p.Results.stimulationIndices;
  stimulationTimeValues = p.Results.stimulationTimeIndices;
  svdFilter = p.Results.svdFilter;
  printTensorStats = p.Results.printStats;
  locCount = p.Results.locCount;
  fwhm = p.Results.fwhm;
  nLocalMax = p.Results.nLocalMax;
  locParams = p.Results.locParams;

  [~, ~, frames] = size(img);
  img = double(img);
  rawData = img;
  if ~isreal(rawData)
    img = abs(img);
    rawData = abs(rawData);
  end

  if ~isempty(svdFilter) && length(svdFilter) > 2 && ndims(img) == 3
    [filteredImg, ~] = filterSVD(img, svdFilter);
    img = sqrt(sum(abs(filteredImg) .^ 2, 3));
  end
  
  if p.Results.depthScatterNorm
    img = depthScatterNormalization(img);
  end

  if p.Results.backNorm
    img = backgroundScatterNormalization(img);
  end

  if locCount > 0
    locParam.fwhm = [fwhm, fwhm];
    locParam.NLocalMax = nLocalMax;
    locParam.InterpMethod = 'spline';
    locParam.LocMethod = 'radial';
    locParam.numberOfParticles = locCount;
    localizedCoords = ULM_localization2D(img, locParam);
    localizedCoordCount = size(localizedCoords, 1);
    if localizedCoordCount == 0
      warning('No Localized Coordinates Found. Localization Failed.');
      localizedCoordCell = {};
    else
      fprintf('Localized %d Points Across %d Frames\n', localizedCoordCount, size(img, 3));
      localizedCoordCell = reshapeLocalizedCoords(localizedCoords);
    end
  else
    localizedCoordCount = 0;
    localizedCoordCell = {};
  end
  
  switch scale
    case 'linear'
    case 'log'
      img = 20 * log10(abs(img));
      img(img < 0) = 0;
    case 'db'
      img = 20 * log10(abs(img));
      img(img == -inf) = 0;
      img = img - max(img(:));
    case 'icg'
      % LINEAR SCALE
    case 'exp'
      img = img .^ expScale;
    otherwise
      error('Unsupported Scale: %s', scale);
  end
  cLimits = [min(img(:)) max(img(:))];

  if ndims(img) == 3 && printTensorStats
    printTensorStatistics(img);
  end

  if ~isempty(p.Results.filtering)
    filterOptions = p.Results.filtering;
    img = filterImageDisplay(img, filterOptions{1}(1), filterOptions{2});
  end

  if ~isempty(p.Results.colorLimits)
    userColorLimits = p.Results.colorLimits;
    cLimits = userColorLimits;
  end

  screenSize = get(0, 'ScreenSize');
  screenWidth = screenSize(3);
  screenHeight = screenSize(4);
  figWidth = 800;
  figHeight = 700;
  figX = (screenWidth - figWidth) / 2;
  figY = (screenHeight - figHeight) / 2;

  fig = uifigure('Name', propertyTitle, 'Color', [0.9412 0.9412 0.9412], ...
                'Position', [figX figY figWidth figHeight]);
  movegui(fig, 'center');
  ax = uiaxes(fig);
  ax.Position = [16 165 768 511];
  ax.FontName = 'Arial';
  ax.FontSize = 14;

  if ~isempty(axesCoords)
    hImg = imagesc(ax, axesCoords{1}, axesCoords{2}, img(:,:,1));
  else
    hImg = imagesc(ax, img(:,:,1));
  end
  axis(ax, 'image');
  colormap(ax, cmap);
  colorbar(ax);
  clim(ax, cLimits);

  if localizedCoordCount > 0
    frameCoords = localizedCoordCell{1};
    hold(ax, 'on');
    scatterHandle = scatter(ax, frameCoords(:, 2), frameCoords(:, 1), 'filled', locParams{1}, 'Marker', 'o');
    scatterHandle.SizeData = fwhm * locParams{2};
    hold(ax, 'off');
  end
  
  if ~isempty(p.Results.xlabel)
    xlabel(ax, p.Results.xlabel);
  end
  if ~isempty(p.Results.zlabel)
    ylabel(ax, p.Results.zlabel);
  end
  if ~isempty(p.Results.axesTitle)
    title(ax, p.Results.axesTitle, 'FontSize', 16, 'FontName', 'Arial', 'Interpreter', 'tex');
  elseif ~isempty(frameTitles)
    title(ax, frameTitles{1}, 'FontSize', 16, 'FontName', 'Arial', 'Interpreter', 'tex');
  end

  if ~isempty(p.Results.mask)
    maskInput = p.Results.mask;
    if isequal(size(maskInput), size(img(:, :, 1)))
      hold(ax, 'on');
      [B, ~] = bwboundaries(maskInput);
      for k = 1:length(B)
         boundary = B{k};
         if ~isempty(axesCoords)
           xVec = axesCoords{1};
           zVec = axesCoords{2};
           plot(ax, xVec(round(boundary(:, 2))), zVec(round(boundary(:, 1))), 'w', 'LineWidth', 2, 'Tag', 'maskBoundary');
         else
           plot(ax, boundary(:, 2), boundary(:, 1), 'w', 'LineWidth', 2, 'Tag', 'maskBoundary');
         end
      end
      hold(ax, 'off');
    else
      warning('Provided mask size does not match image size. Mask overlay skipped.');
    end
  end
  colorLabel = uilabel(fig);
  colorLabel.Text = sprintf('Color: [%.3f, %.3f]', cLimits(1), cLimits(2));
  colorLabel.Position = [16 135 300 22];
  colorLabel.FontWeight = 'bold';
  colorLabel.FontSize = 14;
  
  if ~isempty(p.Results.hardColorLimits)
     sliderRange = p.Results.hardColorLimits;
     ticksValues = round((sliderRange(1):10:sliderRange(2)) / 10) * 10;
     colorSlider = uislider(fig, 'range', 'Orientation', 'horizontal', 'Limits', sliderRange, 'Value', cLimits, ...
                          'MajorTicks', ticksValues); 
  elseif exist('userColorLimits', 'var')
    sliderRange = [min(img(:)) max(img(:))];
    sliderRange(1) = min(sliderRange(1), userColorLimits(1));
    sliderRange(2) = max(sliderRange(2), userColorLimits(2));
    
    ticksValues = round((sliderRange(1):10:sliderRange(2)) / 10) * 10;
    colorSlider = uislider(fig, 'range', 'Orientation', 'horizontal', 'Limits', sliderRange, 'Value', userColorLimits, ...
                          'MajorTicks', ticksValues);
  elseif strcmp(propertyName, 'relativePowerDopplerStack') || strcmp(propertyName, 'relativeRegisteredPowerDopplerStack')
    colorSlider = uislider(fig, 'range', 'Orientation', 'horizontal', 'Limits', [-50 200], 'Value', [-30 120], ...
                          'MajorTicks', [-50 0 50 100 150 200]);
  elseif strcmp(propertyName, 'powerDopplerStack') && strcmp(scale, 'db')
    colorSlider = uislider(fig, 'range', 'Orientation', 'horizontal', 'Limits', [-100 0], 'Value', [-80 0], ...
                          'MajorTicks', [-100 -75 -50 -25 0]);
  else
    colorSlider = uislider(fig, 'range', 'Orientation', 'horizontal', 'Limits', cLimits, 'Value', cLimits, ...
                          'MajorTicks', round(linspace(cLimits(1), cLimits(2), 4)));
  end
  colorSlider.Position = [16 124 768 3];
  colorSlider.ValueChangingFcn = @(s, evt) updateCLim(evt.Value);
  if frames > 1
    frameSlider = uislider(fig, 'Limits', [1 frames], 'Value', 1);
    frameSlider.Position = [220 84 564 3];
    tickCandidates = 10:10:frames;
    if isempty(tickCandidates)
        frameTicks = unique([1, frames]);
    elseif numel(tickCandidates) > 8
        step = ceil(numel(tickCandidates) / 8);
        frameTicks = tickCandidates(1:step:end);
    else
        frameTicks = tickCandidates;
    end
    frameSlider.MajorTicks = frameTicks;
    frameLabel = uilabel(fig);
    frameLabel.Text = 'Frame: 1';
    frameLabel.Position = [16 75 120 22];
    frameLabel.FontWeight = 'bold';
    frameLabel.FontSize = 14;
    
    prevBtn = uibutton(fig, 'push', 'Text', '<', 'Position', [140 75 30 22], ...
                       'ButtonPushedFcn', @(btn, event) changeFrame(-1));
    nextBtn = uibutton(fig, 'push', 'Text', '>', 'Position', [175 75 30 22], ...
                       'ButtonPushedFcn', @(btn, event) changeFrame(1));

    frameSlider.ValueChangedFcn = @(s, e) updateFrame();
  end

  exportFrameBtn = uibutton(fig, 'push', 'Text', 'Export Json', 'Position', [16 10 90 22], ...
                       'ButtonPushedFcn', @(btn, event) exportFigureOptions());
  
  exportRawDataBtn = uibutton(fig, 'push', 'Text', 'Export Data', 'Position', [116 10 90 22], ...
                             'ButtonPushedFcn', @(btn, event) exportRawDataToWorkspace());

  figBtn = uibutton(fig, 'push', 'Text', 'Save Figure', 'Position', [216 10 90 22], ...
                    'ButtonPushedFcn', @(btn, event) saveFigure());

  vidBtn = uibutton(fig, 'push', 'Text', 'Save Video', 'Position', [316 10 90 22], ...
                    'ButtonPushedFcn', @(btn, event) saveVideo());

  svgBtn = uibutton(fig, 'push', 'Text', 'Save PDF', 'Position', [416 10 85 22], ...
                    'ButtonPushedFcn', @(btn, event) saveImage('pdf'));
  
  jpegBtn = uibutton(fig, 'push', 'Text', 'Save JPEG', 'Position', [511 10 85 22], ...
                     'ButtonPushedFcn', @(btn, event) saveImage('jpeg'));

  statsBtn = uibutton(fig, 'push', 'Text', 'Print Statistics', 'Position', [606 10 120 22], ...
                      'ButtonPushedFcn', @(btn, event) printStatistics());

  if ismember(select, {'rectangle', 'polygon', 'mask'})
    selectRegionBtn = uibutton(fig, 'push', 'Text', 'Select Region', 'Position', [16 48 120 22], ...
                              'ButtonPushedFcn', @onSelectRegion);
    selectRegionBtn.BackgroundColor = [0.2 0.6 0.9];
    selectRegionBtn.FontColor = [1 1 1];
    fig.CloseRequestFcn = @(src, ~) closeAndResume(src);
    uiwait(fig);
  end

  roiCell = {select, roi, tag};

  function closeAndResume(src)
    uiresume(src);
    delete(src);
  end

  function onSelectRegion(~, ~)
    selectRegionBtn.Enable = 'off';
    selectRegionBtn.Text = 'Drawing...';
    drawnow;
    switch select
      case 'rectangle'
        roi = zeros(1, 4);
        rect = drawrectangle(ax);
        wait(rect);
        if isvalid(rect)
          roi(1) = round(rect.Position(2));
          roi(2) = round(rect.Position(2)) + round(rect.Position(4)) - 1;
          roi(3) = round(rect.Position(1));
          roi(4) = round(rect.Position(1)) + round(rect.Position(3)) - 1;
          fprintf('Selected Rectangle ROI - %s\n', tag)
          if seriesPlot
            [timeSeries, averagedSeries, seriesStd] = rectangleSeriesGeneration(rawData, roi);
            plotSeriesData(timeSeries, averagedSeries, seriesStd);
          end
        end
        clear rect;
      case 'polygon'
        poly = drawpolygon(ax);
        wait(poly);
        if isvalid(poly)
          mask = createMask(poly, hImg);
          roi = mask;
          fprintf('Selected Polygon ROI - %s\n', tag)
          if seriesPlot
            [timeSeries, averagedSeries, seriesStd] = polygonSeriesGeneration(rawData, roi);
            plotSeriesData(timeSeries, averagedSeries, seriesStd);
          end
        end
        clear poly;
      case 'mask'
        impoly = drawfreehand(ax);
        wait(impoly);
        if isvalid(impoly)
          mask = createMask(impoly, hImg);
          roi = mask;
          fprintf('Selected Mask ROI - %s\n', tag)
          if seriesPlot
            warning('Mask ROI Series Plotting Not Implemented Yet.');
          end
        end
        clear impoly;
    end
    selectRegionBtn.Enable = 'on';
    selectRegionBtn.Text = 'Select Region';
    uiresume(fig);
  end

  function changeFrame(delta)
    if frames > 1
      newVal = frameSlider.Value + delta;
      newVal = max(1, min(frames, newVal));
      frameSlider.Value = newVal;
      updateFrame();
    end
  end

  function printStatistics()
    if frames > 1
      idx = round(frameSlider.Value);
    else
      idx = 1;
    end
    currentFrame = rawData(:,:,idx);
    meanVal = mean(currentFrame(:));
    medianVal = median(currentFrame(:));
    stdVal = std(currentFrame(:));
    
    fprintf('Statistics for Frame %d:\n', idx);
    fprintf('  Mean: %.4f\n', meanVal);
    fprintf('  Median: %.4f\n', medianVal);
    fprintf('  Standard Deviation: %.4f\n', stdVal);
  end

  function updateFrame()
    idx = round(frameSlider.Value);
    hImg.CData = img(:,:,idx);
    frameLabel.Text = sprintf('Frame: %d', idx);
    if ~isempty(stimulationIndices) && length(stimulationIndices) >= idx
      if stimulationIndices(idx)
        title(ax, sprintf('Block: [%d/%d] Time: %.2fs CO2 ON', idx, size(stimulationIndices, 1), stimulationTimeValues(idx)), 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0 0.6 0], 'FontName', 'Arial');
      else
        title(ax, sprintf('Block: [%d/%d] Time: %.2fs CO2 OFF', idx, size(stimulationIndices, 1), stimulationTimeValues(idx)), 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'r', 'FontName', 'Arial');
      end
    elseif ~isempty(frameTitles) && length(frameTitles) >= idx
      title(ax, frameTitles{idx}, 'FontSize', 20, 'FontName', 'Arial', 'Interpreter', 'tex', 'Color', 'k', 'FontWeight', 'normal');
    end
    
    if localizedCoordCount > 0
      if exist('scatterHandle', 'var') && isvalid(scatterHandle)
        delete(scatterHandle);
      end
      frameCoords = localizedCoordCell{idx};
      hold(ax, 'on');
      scatterHandle = scatter(ax, frameCoords(:, 2), frameCoords(:, 1), 'filled', locParams{1}, 'Marker', 'o');
      scatterHandle.SizeData = fwhm * locParams{2};
      hold(ax, 'off');
    end
    
    hBoundaries = findobj(ax, 'Tag', 'maskBoundary');
    if ~isempty(hBoundaries)
      uistack(hBoundaries, 'top');
    end
  end

  function updateCLim(val)
    caxis(ax, val);
    colorLabel.Text = sprintf('Color: [%.3f, %.3f]', val(1), val(2));
  end

  function saveImage(format)
    filter = {'*.*', 'All Files'};
    if strcmp(format, 'pdf')
      filter = {'*.pdf', 'PDF Files (*.pdf)'};
      defaultName = 'image.pdf';
    elseif strcmp(format, 'jpeg')
      filter = {'*.jpg;*.jpeg', 'JPEG Files (*.jpg, *.jpeg)'};
      defaultName = 'image.jpg';
    end
    [file, path] = uiputfile(filter, 'Save Image As', defaultName);
    if isequal(file, 0) || isequal(path, 0)
      return;
    end
    savePath = fullfile(path, file);
    try
      if strcmp(format, 'pdf')
        exportgraphics(fig, savePath, 'ContentType', 'vector');
      else
        exportgraphics(fig, savePath, 'Resolution', 300);
      end
      fprintf('Image Saved to: %s\n', savePath);
      uialert(fig, sprintf('Image successfully saved to: %s', file), 'Export Success', 'Icon', 'success');
    catch ME
      uialert(fig, sprintf('Failed to Save Image: %s', ME.message), 'Export Note', 'Icon', 'info');
    end
  end

  function saveVideo()
    [file, path] = uiputfile({'*.mp4', 'MP4 Files (*.mp4)'; '*.avi', 'AVI Files (*.avi)'}, 'Save Video As');
    if isequal(file, 0) || isequal(path, 0)
      return;
    end
    savePath = fullfile(path, file);
    [~, ~, ext] = fileparts(file);
    profiles = VideoWriter.getProfiles();
    profileNames = {profiles.Name};
    if strcmpi(ext, '.mp4') && ismember('MPEG-4', profileNames)
      vid = VideoWriter(savePath, 'MPEG-4');
    elseif strcmpi(ext, '.mp4')
      [pathstr, name] = fileparts(savePath);
      savePath = fullfile(pathstr, [name '.avi']);
      fprintf('MP4 is not supported on this system. Saving as AVI (Motion JPEG) instead.\n');
      vid = VideoWriter(savePath, 'Motion JPEG AVI');
    elseif ismember('Motion JPEG AVI', profileNames)
      vid = VideoWriter(savePath, 'Motion JPEG AVI');
    else
      vid = VideoWriter(savePath, 'Uncompressed AVI');
    end
    vid.FrameRate = 5;
    open(vid);
    currentCLim = colorSlider.Value;
    for frameIdx = 1:size(img, 3)
      if frames > 1
        frameSlider.Value = frameIdx;
        frameLabel.Text = sprintf('Frame: %d', frameIdx);
      end
      hImg.CData = img(:,:,frameIdx);
      clim(ax, currentCLim);
      if ~isempty(frameTitles) 
        title(ax, frameTitles{frameIdx}, 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Arial', 'Interpreter', 'tex');
      elseif ~isempty(stimulationIndices) && length(stimulationIndices) >= frameIdx
        if stimulationIndices(frameIdx)
          title(ax, sprintf('Block [%d/%d] CO2 ON', frameIdx, size(stimulationIndices, 1)), 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0 0.6 0], 'FontName', 'Arial');
        else
          title(ax, sprintf('Block [%d/%d] CO2 OFF', frameIdx, size(stimulationIndices, 1)), 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'r', 'FontName', 'Arial');
        end
      else
        title(ax, sprintf('Frame Index: [%d/%d]', frameIdx, size(img, 3)), 'FontSize', 20, 'FontWeight', 'bold', 'FontName', 'Arial', 'Interpreter', 'tex');
      end
      hBoundaries = findobj(ax, 'Tag', 'maskBoundary');
      if ~isempty(hBoundaries)
        uistack(hBoundaries, 'top');
      end
      drawnow;
      frame = getframe(fig);
      writeVideo(vid, frame);
    end
    close(vid);
    updateFrame();
  end

  % DEPRECATED FUNCTION
  function exportFrameToWorkspace()
    prompt = {'Enter variable name:'};
    dlgtitle = 'Export Current Frame to Workspace';
    dims = [1 35];
    definput = {'exportedRawImg'};
    answer = inputdlg(prompt, dlgtitle, dims, definput);
    if ~isempty(answer)
      varName = answer{1};
      if isvarname(varName)
        if frames > 1
          idx = round(frameSlider.Value);
          exportData = img(:,:,idx);
        else
          exportData = img;
        end
        assignin('base', varName, exportData);
        fprintf('Exported Current Frame to Workspace - "%s"\n', varName);
        uialert(fig, sprintf('Frame exported as "%s"', varName), 'Export Success', 'Icon', 'success');
      else
        uialert(fig, 'Invalid variable name.', 'Export Note', 'Icon', 'info');
      end
    end
  end

  function saveFigureOptions()
    [file, path] = uiputfile('*.json', 'Save Figure Options - JSON');
    if isequal(file, 0) || isequal(path, 0)
      return;
    end
    savePath = fullfile(path, file);
    try
      options = p.Results;
      if exist('colorSlider', 'var') && isvalid(colorSlider)
        options.colorLimits = colorSlider.Value;
      end
      fieldsToRemove = {'mask', 'axesCoords'}; 
      for i = 1:length(fieldsToRemove)
        if isfield(options, fieldsToRemove{i})
          options = rmfield(options, fieldsToRemove{i});
        end
      end
      jsonStr = jsonencode(options);
      fid = fopen(savePath, 'w');
      if fid == -1
        error('Cannot Create File: %s', savePath);
      end
      fprintf(fid, '%s', jsonStr);
      fclose(fid);
      uialert(fig, sprintf('Options saved to: %s', file), 'Export Success', 'Icon', 'success');
    catch ME
      uialert(fig, sprintf('Failed to Save Figure Options: %s', ME.message), 'Export Note', 'Icon', 'info');
    end
  end

  function exportRawDataToWorkspace()
    prompt = {'Enter variable name:'};
    dlgtitle = 'Export Raw Data to Workspace';
    dims = [1 35];
    definput = {'exportedRawImg'};
    answer = inputdlg(prompt, dlgtitle, dims, definput);
    if ~isempty(answer)
      varName = answer{1};
      if isvarname(varName)
        assignin('base', varName, rawData);
        fprintf('Exported Raw Data to Workspace - "%s"\n', varName);
        uialert(fig, sprintf('Data stack exported as "%s"', varName), 'Export Success', 'Icon', 'success');
      else
        uialert(fig, 'Invalid variable name.', 'Export Note', 'Icon', 'info');
      end
    end
  end

  function saveFigure()
    [file, path] = uiputfile('*.fig', 'Save Figure As');
    if isequal(file, 0) || isequal(path, 0)
      return;
    end
    savePath = fullfile(path, file);
    try
      saveas(fig, savePath);
      fprintf('Figure Saved to: %s\n', savePath);
      uialert(fig, sprintf('Figure successfully saved to: %s', file), 'Export Success', 'Icon', 'success');
    catch ME
      uialert(fig, sprintf('Failed to Save Figure: %s', ME.message), 'Export Note', 'Icon', 'info');
    end
  end

end

function plotSeriesData(timeSeries, averagedSeries, seriesStd, varargin)
  fig = figure;
  set(fig, 'Position', [400, 400, 1000, 350]);
  movegui(fig, 'center');
  t = 1:numel(timeSeries);
  t = t(:)'; 
  timeSeries = timeSeries(:)'; 
  averagedSeries = averagedSeries(:)'; 
  seriesStd = seriesStd(:)';
  plot(t, timeSeries, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5, 'DisplayName', 'Raw Time Series');
  hold on;
  plot(t, averagedSeries, 'b', 'LineWidth', 5, 'DisplayName', 'Processed Time Series');
  yLower = averagedSeries - seriesStd;
  yUpper = averagedSeries + seriesStd;
  fill([t, fliplr(t)], [yUpper, fliplr(yLower)], 'b', ...
       'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
  hold off;
  xlabel('Frame Index', 'FontWeight', 'bold', 'FontSize', 14);
  ylabel('Intensity', 'FontWeight', 'bold', 'FontSize', 14);
  titleStr = sprintf('Mean - %.1f, Std - %.1f, Range - [%.1f, %.1f], COV - %.1f%%', ...
      mean(timeSeries), std(timeSeries), min(timeSeries), max(timeSeries), ...
      std(timeSeries) / mean(timeSeries) * 100);
  title(titleStr, 'FontWeight', 'bold', 'FontSize', 16);
  legend('show');
  grid on;
end

function localStimFill(timeSeries, stimSections, stimStruct, y1, y2)
  for k = stimSections
    r = stimStruct{k};
    r(1) = max(1, r(1)); r(2) = min(numel(timeSeries), r(2));
    x1 = timeSeries(r(1)); x2 = timeSeries(r(2));
    fill([x1 x2 x2 x1], [y1 y1 y2 y2], [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
  end
end

function normalizedData = depthScatterNormalization(rawData)
  if ndims(rawData) > 2
    normalizedData = zeros(size(rawData));
    for idx = 1:size(rawData, 3)
      normalizedData(:,:,idx) = depthScatterNormalization2D(rawData(:,:,idx));
    end
  elseif ndims(rawData) == 2
    normalizedData = depthScatterNormalization2D(rawData);
  end
end

function normalizedData = depthScatterNormalization2D(rawData)
  rowNormFactor = mean(rawData(:, 1:5), 2);
  normalizedData = rawData ./ rowNormFactor * rowNormFactor(1);
end

function normalizedData = backgroundScatterNormalization(rawData)
  if ndims(rawData) == 3
    backgroundFactor = mean(rawData(1:5, end-5:end, :), [1, 2]); 
  else
    backgroundFactor = mean(rawData(1:5, end-5:end), [1, 2]);
  end
  normalizedData = rawData ./ backgroundFactor;
end

function printTensorStatistics(dataTensor)
  fprintf('Frame Count: %d\n', size(dataTensor, 3));
  fprintf('Mean: %f\n', mean(dataTensor, 'all'));
  fprintf('Median: %f\n', median(dataTensor, 'all'));
  fprintf('Std: %f\n', std(dataTensor, 0, 'all'));
  fprintf('Range: [%f, %f]\n', min(dataTensor, [], 'all'), max(dataTensor, [], 'all'));
end
