function [stat, info] = column_stat(x_px, y_value, which_stat, min_frame)
%COLUMN_STAT  One number per one-pixel column of the x series.
%   Every thickness in this pipeline is a difference of integer row indices, so a
%   one-pixel column of x holds exactly one attainable value and grouping by it is
%   exact rather than a binning choice.
%
% IN   x_px       T x 1 double px  the column variable, integer valued
%      y_value    T x 1 double     what is summarised in each column, any unit
%      which_stat 1 x 1 str        "median" | "centroid" | "mode"
%      min_frame  1 x 1 int        frames a column needs before it gives a point
% OUT  stat  1 x nc double  the statistic, in y_value's unit
%       info  struct
%         .x_px     1 x nc double px  the column each entry belongs to
%         .n_frame  1 x nc double     frames behind each entry
%         .n_column 1 x 1 int         columns that cleared min_frame
%
%   "median" IS READ BETWEEN THE LATTICE POINTS. y is usually also a multiple of
%   the pixel size, so a plain median lands on the lattice and a difference of two
%   of them can come out as exactly zero for no reason but the grid. The tie block
%   is interpolated instead. see CLAUDE_LOG.md
    arguments
        x_px       (:,1) double
        y_value    (:,1) double
        which_stat (1,1) string {mustBeMember(which_stat, ["median", "centroid", "mode"])}
        min_frame  (1,1) double {mustBeInteger, mustBePositive}
    end

usable = isfinite(x_px) & isfinite(y_value);
column_value = unique(x_px(usable));
stat = nan(1, numel(column_value));
n_frame = zeros(1, numel(column_value));
for k = 1:numel(column_value)
    inside = usable & x_px == column_value(k);
    n_frame(k) = sum(inside);
    if n_frame(k) < min_frame
        continue
    end
    y_here = y_value(inside);
    switch which_stat
        case "centroid"
            stat(k) = mean(y_here);
        case "mode"
            stat(k) = mode(y_here);
        otherwise
            stat(k) = median_between(y_here);
    end
end

keep = isfinite(stat);
stat = stat(keep);
info.x_px = column_value(keep)';
info.n_frame = n_frame(keep);
info.n_column = sum(keep);
end

function value = median_between(y_here)
    % the median, then spread across whatever tie block it landed in, so the answer
    % is not pinned to the lattice the integer boundaries allow
    y_here = sort(y_here);
    n_total = numel(y_here);
    value = median(y_here);
    step = min(diff(unique(y_here)));
    if isempty(step) || step <= 0
        return
    end
    tied = abs(y_here - value) < step / 2;
    n_tied = sum(tied);
    if n_tied < 2
        return
    end
    n_below = sum(y_here < value - step / 2);
    value = value - step / 2 + step * (n_total / 2 - n_below) / n_tied;
end
