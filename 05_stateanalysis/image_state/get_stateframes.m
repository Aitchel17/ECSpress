function transition_frames = get_stateframes(state_frameidx, frames, pad)
%GET_STATEFRAMES  Cut each state bout out of the image stack, with optional margin.
%
%   transition_frames = get_stateframes(state_frameidx, frames)
%   transition_frames = get_stateframes(state_frameidx, frames, pad)
%
%   state_frameidx : B x 2 matrix, each row [start end] frame index of one bout.
%   frames         : H x W x N image stack.
%   pad            : frames to extend each bout by. Scalar -> same before/after;
%                    [pre post] -> asymmetric. Default 0. Clamped to [1, N] so
%                    bouts near the recording ends are simply not extended past them.
%
%   Returns a 1 x B cell; each cell is H x W x (bout length + padding).

    arguments
        state_frameidx double
        frames         {mustBeNumeric}
        pad (1,:) double = 0
    end

    if isscalar(pad)
        pre = pad;  post = pad;
    else
        pre = pad(1);  post = pad(2);
    end

    nframes = size(frames, 3);

    transition_frames = cell(1, size(state_frameidx, 1));
    for bout_idx = 1:size(state_frameidx, 1)
        fidx = state_frameidx(bout_idx, :);
        f1 = max(1,       fidx(1) - pre);
        f2 = min(nframes, fidx(2) + post);
        transition_frames{bout_idx} = frames(:, :, f1:f2);
    end
end
