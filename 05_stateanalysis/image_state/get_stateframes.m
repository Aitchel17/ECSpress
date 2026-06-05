function transition_frames = get_stateframes(state_frameidx,frames)
    %  get cell of frames
    
    transition_frames = {};
    for bout_idx = 1:size(state_frameidx,1)
            fidx = state_frameidx(bout_idx,:);
            transition_frames{bout_idx} = frames(:,:,fidx(1):fidx(2));
    end

end

