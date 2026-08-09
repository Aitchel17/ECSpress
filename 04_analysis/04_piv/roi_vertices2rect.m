function roirect = roi_vertices2rect(vertices)
%ROI_VERTICES2RECT  Convert ROI vertices to a bounding-box rect.
%   Takes the vertex list produced by roi_handle / roiSelector and returns the
%   [xmin ymin width height] rect expected by PIVlab's piv_preprocess,
%   rounded to whole pixels.
%
% IN   vertices  Nx2 [x y] ROI vertices; empty in gives empty out
% OUT  roirect   1x4 [xmin ymin width height] in integer pixels

    % 1. Empty ROI passes through
    if isempty(vertices)
        roirect = [];
        return;
    end

    % 2. Bounding box of the vertices
    xmin = min(vertices(:, 1));
    ymin = min(vertices(:, 2));
    xmax = max(vertices(:, 1));
    ymax = max(vertices(:, 2));

    width  = xmax - xmin;
    height = ymax - ymin;

    % 3. Round to pixel coordinates
    roirect = round([xmin, ymin, width, height]);
end
