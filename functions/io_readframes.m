function zstack = io_readframes(mobj, imgch, zrange, batchSize, verbose)
    % Function to read frames from ActiveX-based .mdf file in batches
    % Inputs:
    %   mobj - ActiveX object with ReadFrame method
    %   imgch - Image channel to read
    %   zrange - Range of frames [start, end]
    %   batchSize - (Optional) Number of frames per batch (default: 1000)
    %   verbose - (Optional) Boolean to enable/disable progress messages
    % Output:
    %   zstack - 3D stack of frames

    if nargin < 4 || isempty(batchSize)
        batchSize = 1000; % Default batch size
    end
    if nargin < 5
        verbose = true; % Default to verbose output
    end
    disp(['Loading Image ch',string(imgch)])

    % Determine the number of frames and batches
    totalFrames = zrange(2) - zrange(1) + 1;
    numBatches = ceil(totalFrames / batchSize);

    % Pre-allocate zstack as a cell array to hold batch results
    zstack = cell(1, numBatches);

    % Process each batch
    tic;
    for b = 1:numBatches
        % Compute frame range for the current batch
        batchStart = zrange(1) + (b - 1) * batchSize;
        batchEnd = min(zrange(2), batchStart + batchSize - 1);

        if verbose
            fprintf('Reading frames batch %d/%d (frames %d to %d)...\n', ...
                b, numBatches, batchStart, batchEnd);
        end

        % Read the batch of frames
        zstack{b} = io_readframes_simple(mobj, imgch, [batchStart, batchEnd], false);
    end
    toc;

    % Combine batches into a single 3D stack
    if verbose
        fprintf('Combining batches into final Z-stack...\n');
    end
    zstack = cat(3, zstack{:});


end

function zstack = io_readframes_simple(mobj, imgch, zrange, verbose)
    % Helper function to sequentially read frames
    % Inputs:
    %   mobj - ActiveX object
    %   imgch - Image channel
    %   zrange - Range of frames [start, end]
    %   verbose - Enable/disable logging
    % Outputs:
    %   zstack - 3D image stack

    start_z = zrange(1);
    end_z = zrange(2);
    lenz = end_z - start_z + 1;

    % Read first frame to get dimensions
    sampleFrame = mobj.ReadFrame(imgch, start_z)';
    [height, width] = size(sampleFrame);
    zstack = zeros(height, width, lenz, 'like', sampleFrame);

    % Read each frame sequentially
    for idx = 1:lenz
        z = start_z + idx - 1; % Compute absolute frame index
        if verbose
            fprintf('Reading frame %d of %d...\n', idx, lenz);
        end
        zstack(:, :, idx) = mobj.ReadFrame(imgch, z)';
    end
end
