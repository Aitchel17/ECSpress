function channelData = analyze_readtiff(folderdirectory,namingpattern)
tic
% find .tif file
    channelDir = dir(fullfile(folderdirectory, namingpattern));
    if length(channelDir) ~= 1
        fprintf('%s file not exist or plural num: %d\n', namingpattern, length(channelDir));
        channelData = [];
        return;
    end
   
    filePath = fullfile(folderdirectory, channelDir.name);
    
% load metadata info    
    info = imfinfo(filePath);
    numFrames = numel(info);

    % Determine data type
    dataType = info(1).BitDepth;
    if dataType == 16
        dataClass = 'uint16';
    elseif dataType == 8
        dataClass = 'uint8';
    else
        error('Unsupported BitDepth: %d', dataType);
    end

% Loading 
    % Preallocate the array
    channelData = zeros(info(1).Height, info(1).Width, numFrames, dataClass);

    % Progress bar update rule
    updateInterval = max(1, round(numFrames / 100)); 
    D = parallel.pool.DataQueue;
    progress = 0;
    h = waitbar(0, sprintf('Loading %s...', namingpattern));
    afterEach(D, @(~) updateWaitbar());

    % parallel loading
    parfor idx = 1:numFrames
        tempTiff = Tiff(filePath, 'r');
        cleanup = onCleanup(@() tempTiff.close());         % object that trigger .close() when it destroied
        tempTiff.setDirectory(idx);
        channelData(:, :, idx) = tempTiff.read();
        if mod(idx, updateInterval) == 0
            send(D, idx);
        end
    end
    close(h);

    function updateWaitbar()
        progress = progress + updateInterval;
        waitbar(min(progress / numFrames, 1), h);
    end
toc
end
