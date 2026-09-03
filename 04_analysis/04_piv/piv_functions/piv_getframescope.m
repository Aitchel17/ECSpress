function [stack_span, span] = piv_getframescope(stack, from, to, halfwin)
%PIV_GETFRAMESCOPE  The frames one event needs, cut out of the recording.
%   from-halfwin : to+halfwin. An endpoint too close to either end is refused rather
%   than moved, because moving it would change which frames the event is.
%
% IN   stack       H x W x T numeric   the recording, before any PIV preprocessing
%      from        1 x 1 int           the earlier endpoint, a frame index into stack
%      to          1 x 1 int           the later one
%      halfwin     1 x 1 int           frames either side of each endpoint
% OUT  stack_span  H x W x n_span      stack(:, :, span)
%      span        1 x n_span int      the frame indices cut, for the record
    arguments
        stack   {mustBeNumeric, mustBeNonempty}
        from    (1,1) double {mustBeInteger, mustBePositive}
        to      (1,1) double {mustBeInteger, mustBePositive}
        halfwin (1,1) double {mustBeInteger, mustBeNonnegative}
    end
    n_frame = size(stack, 3);
    span = from - halfwin : to + halfwin;
    if span(1) < 1 || span(end) > n_frame
        error('piv_getframescope:margin', ...
            'event %d-%d needs frames %d-%d for halfwin %d, and the recording is %d long', ...
            from, to, span(1), span(end), halfwin, n_frame);
    end
    stack_span = stack(:, :, span);
end
