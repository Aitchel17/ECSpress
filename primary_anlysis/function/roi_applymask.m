function maskedstack = roi_applymask(stack, mask)
    [rows, cols] = find(mask);
    ymin = min(rows);    ymax = max(rows);
    xmin = min(cols);    xmax = max(cols);
    %
    cropstack = stack(ymin:ymax, xmin:xmax, :);
    cropmask = mask(ymin:ymax, xmin:xmax);
    maskedstack = double(cropstack);
    maskedstack(~repmat(cropmask, 1, 1, size(cropstack,3))) = NaN;
end
