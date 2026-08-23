function [code, info] = state_perframe(t_axis, behav_state, binwidth_sec, edge_frame)
%STATE_PERFRAME  The arousal score laid on the imaging frame axis.
%   The scoring lives on its own coarse grid -- binwidth_sec seconds per bin --
%   and every frame falls inside one of those bins. This says which.
%
%   IT RETURNS THE CODE, NOT A NAME. statecodes is the scoring's own lookup and
%   differs between scoring generations, so the caller compares against the codes
%   it read from the same file rather than trusting a name written here.
%
% IN   t_axis        1 x T float s    frame times, from the recording
%      behav_state   1 x B double     the score, one entry per bin
%      binwidth_sec  1 x 1 float s    seconds per scoring bin
%      edge_frame    1 x 1 int        frames either side of a label change to flag.
%                                     0 flags only the change frames themselves
% OUT  code  T x 1 double  the scoring code for each frame
%      info  struct
%        .near_change   T x 1 bool  true = within edge_frame of a change
%        .n_change      1 x 1 int   label changes in the record
%        .frac_flagged  1 x 1 float what fraction near_change covers
%
%   WHY THE BAND IS A FLAG AND NOT A CUT. Mislabelling near a change makes two
%   state cells more alike, so it shrinks a difference and cannot manufacture one.
%   Dropping the band therefore protects against a false positive at the cost of
%   sessions, and whether that trade is worth taking depends on what is being
%   claimed -- which is the caller's business. see FINDINGS.md
    arguments
        t_axis       (1,:) double
        behav_state  (1,:) double
        binwidth_sec (1,1) double {mustBePositive}
        edge_frame   (1,1) double {mustBeInteger, mustBeNonnegative}
    end

% a frame lands in the bin its own time falls in; the last bin absorbs any frame
% past the end of the scoring, which happens when the recording outruns it
bin_index = floor(t_axis(:) / binwidth_sec) + 1;
bin_index = min(bin_index, numel(behav_state));
code = behav_state(bin_index)';

changed = [true; diff(code) ~= 0];
window = 2 * edge_frame + 1;
info.near_change = movmax(double(changed), window) > 0;
info.n_change = sum(changed) - 1;
info.frac_flagged = mean(info.near_change);
end
