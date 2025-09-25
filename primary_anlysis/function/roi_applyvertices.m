function roi_stack = roi_applyvertices(stack, vertices)
    stack = double(stack);
    xmin = round(min(vertices(:, 1)));
    xmax = round(max(vertices(:, 1)));
    ymin = round(min(vertices(:, 2)));
    ymax = round(max(vertices(:, 2)));

    [cols, rows] = meshgrid(1:size(stack, 2), 1:size(stack, 1));
    mask = inpolygon(cols, rows, vertices(:, 1), vertices(:, 2));

    roi_stack = stack .* mask;
    roi_stack(~mask) = NaN;
    roi_stack = roi_stack(ymin:ymax, xmin:xmax, :);
end