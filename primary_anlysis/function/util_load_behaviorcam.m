function behavior_cam = util_load_behaviordata(folderpath)
%UTIL_LOAD_BEHAVIORDATA Summary of this function goes here
    % load eye.avi, and whisker.avi
    % at the point of 2025.08.17 MScan behavior cam fps has a bug - it does
    % not match with imaging speed 
%   Detailed explanation goes here
    disp ('load behavior data start')
    eye_path = fullfile(folderpath,'eye.avi');
    whisker_path = fullfile(folderpath,'whisker.avi');
    

    % Initialize outputs
    behavior_cam.eye     = NaN;

    if isfile(eye_path)
    else
        disp('eye.avi does not exist')
    end
    
    if isfile(whisker_path)
    else
        disp('whisker.avi does not exist')
    end
    behavior_cam.eye = read_avi(eye_path);
    behavior_cam.whisker = read_avi(whisker_path);
    
    disp ('load behavior data end')

end

function [frames, fps] = read_avi(fpath)
%READ_AVI One-pass AVI loader with preallocation (Duration*FPS estimate).
%   Returns grayscale frames (HxWxN, uint8) and fps.
%   Shows progress in 5% steps.

    [~, tag, ~] = fileparts(fpath);
    frames = [];
    fps    = NaN;

    if ~isfile(fpath)
        fprintf('%s.avi does not exist\n', tag);
        return;
    end

    try
        vr = VideoReader(fpath);
        fps = vr.FrameRate;

        if ~hasFrame(vr)
            warning('%s.avi has no frames.', tag);
            return;
        end

        % --- 첫 프레임 읽고 크기 확인
        firstFrame = readFrame(vr);
        if ndims(firstFrame) == 3
            firstFrame = firstFrame(:,:,1); % grayscale 강제
        end
        [H,W] = size(firstFrame);

        % --- 추정 프레임 수
        nEst = max(1, floor(vr.Duration * vr.FrameRate));

        % --- 배열 preallocate
        frames = zeros(H,W,nEst,'uint8');
        frames(:,:,1) = firstFrame;

        % --- 진행도 설정
        pctNext = 5;

        % --- 읽기 루프
        k = 1;
        while hasFrame(vr)
            k = k + 1;
            f = readFrame(vr);
            if ndims(f) == 3, f = f(:,:,1); end
            if ~isa(f,'uint8'), f = im2uint8(f); end
            if k > size(frames,3)   % overshoot 대비 확장
                frames(:,:,end+1) = 0;
            end
            frames(:,:,k) = f;

            % 진행도 (5% 단위)
            pct = floor((k / nEst) * 100);
            if pct >= pctNext
                fprintf('\r[%s] Progress: %3d%%', tag, pct);
                pctNext = pctNext + 5;
            end
        end

        % --- trim
        frames = frames(:,:,1:k);

        fprintf('\r[%s] Done. %d frames @ %.3f fps\n', tag, k, fps);

    catch ME
        warning('Failed to read %s: %s', tag, ME.message);
        frames = [];
        fps    = NaN;
    end
end

