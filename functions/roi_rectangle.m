function [vertices] = roi_rectangle(stack)
    % Display the first frame of the stack so the user can draw a box on a steady region
    figureHandle = figure(77);
    sliceViewer(stack);
    title('Draw a box around the region to be processed and click Confirm');
    
    % Draw rectangle and get vertices
    thebox = drawrectangle(gca);
    vertices = round(thebox.Vertices);
    
    % Add a Confirm button
    uicontrol('Style', 'pushbutton', ...
              'String', 'Confirm', ...
              'Position', [20, 20, 100, 40], ...
              'Callback', @(src, event) uiresume(gcbf)); % Resume execution when clicked
    
    % Wait for confirmation
    uiwait(figureHandle);
    
    % Close the figure after confirmation
    close(figureHandle);
end
