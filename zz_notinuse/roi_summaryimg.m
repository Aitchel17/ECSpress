function fig = roi_summaryimg(image, vertices, roimod)
    % ROI_SUMMARYIMG Overlays a single ROI on the input image using given vertices.
    %
    % INPUTS:
    %   image    - The 2D image on which the ROI is drawn.
    %   vertices - Nx2 matrix defining the ROI vertices.
    %
    % OUTPUT:
    %   overlay - Image with ROI overlaid.

    % Plot the image and overlay the ROI
    fig = figure();
    im_handle = imshow(image, []);
    hold on;
    
    if strcmp(roimod, 'line')
        plot([vertices(1, 1); vertices(2, 1)], [vertices(1, 2); vertices(2, 2)], 'r-', 'LineWidth', vertices(3, 1)); % Close the polygon
    else
        plot([vertices(:, 1); vertices(1, 1)], [vertices(:, 2); vertices(1, 2)], 'r-', 'LineWidth', 2); % Close the polygon
    end
    hold off;

    % Output the original image (no modifications made to the image itself)
    imcontrast(im_handle)
    overlay = image;
end
