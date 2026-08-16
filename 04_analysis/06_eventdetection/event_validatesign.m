function [n_ok, n_total, ok] = event_validatesign(pol, value)
%EVENT_VALIDATESIGN  Does the measured sign follow the prescribed polarity?  (S1)
%   Counts, and returns the count. It draws nothing.
%
% IN   pol      1 x N string or cellstr, "dilation" / "constriction" / anything
%               else, which is treated as having no prescribed direction
%      value    1 x N float, the measured quantity. Only its SIGN is read
% OUT  n_ok     1 x 1 int, rows whose sign matched
%      n_total  1 x 1 int, rows that HAVE a prescribed direction
%      ok       1 x N logical, which rows matched. A row with no prescribed
%               direction is false here and absent from n_total, so `ok` alone
%               cannot be summed against numel(pol)
%
%   see CLAUDE_LOG.md   why the denominator is not numel(pol), and what the
%                       returned pair is a test statistic against

arguments
    pol   (1,:)
    value (1,:) double
end

pol = string(pol);
if numel(pol) ~= numel(value)
    error('event_validatesign:sizeMismatch', ...
        'pol has %d entries and value has %d.', numel(pol), numel(value));
end

is_dilation     = pol == "dilation";
is_constriction = pol == "constriction";

ok      = (is_dilation & value > 0) | (is_constriction & value < 0);
n_total = nnz(is_dilation | is_constriction);
n_ok    = nnz(ok);
end
