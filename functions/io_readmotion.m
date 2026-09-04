function [disp_row, disp_col, regerror, drift_fps] = io_readmotion(motion_path)
%IO_READMOTION  One *_motion.txt as tissue displacement, sign corrected.
%   Caller: centralize_drift
%   The extractor writes seven lines: a header, driftestimation_fps, a table header,
%   then regerror, diffphase, rowshift and colshift as mat2str of one row each.
%
% IN   motion_path  1 x 1 str    the *_motion.txt
% OUT  disp_row     N x 1 double px, + = tissue moved toward increasing ROW
%      disp_col     N x 1 double px, + = tissue moved toward increasing COLUMN
%      regerror     N x 1 double dimensionless registration residual
%      drift_fps    1 x 1 double rate the drift was estimated at, Hz
%
%   THE SIGN IS FLIPPED HERE, ONCE. The file stores what dft_registration handed to
%   imtranslate, which is the correction that pulls a frame BACK onto the reference,
%   so the tissue moved the other way. PIV reports tissue displacement, and two
%   quantities that point opposite ways in one figure is the error this prevents.
%
%   The file used to label these yerror / xerror / ymotion / xmotion, which named an
%   axis on the two that have none and said y,x for what are row,col. Every file on
%   disk was relabelled on 2026-09-01; a file older than that will not parse.
%   'diffphase' is angle(CCmax), zero for real images, and is not returned.
    text = string(fileread(motion_path));
    line_list = splitlines(text);

    drift_fps = read_scalar(line_list, "driftestimation_fps");
    regerror = read_series(line_list, "regerror");
    disp_row = -read_series(line_list, "rowshift");
    disp_col = -read_series(line_list, "colshift");

    length_list = [numel(regerror), numel(disp_row), numel(disp_col)];
    if numel(unique(length_list)) > 1
        error('io_readmotion:raggedSeries', ...
            'series lengths %s disagree in %s', mat2str(length_list), motion_path);
    end
end

function value = read_scalar(line_list, label)
    hit = find(startsWith(line_list, label + ":"), 1);
    if isempty(hit)
        error('io_readmotion:noLabel', '%s is not in the file', label);
    end
    value = str2double(extractAfter(line_list(hit), label + ":"));
end

function value = read_series(line_list, label)
    hit = find(startsWith(line_list, label + ":"), 1);
    if isempty(hit)
        error('io_readmotion:noLabel', '%s is not in the file', label);
    end
    body = extractAfter(line_list(hit), label + ":");
    body = erase(body, ["[", "]"]);
    value = sscanf(body, '%f');
end
