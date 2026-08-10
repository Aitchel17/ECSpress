function time_idx = statebin2frame(varargin)
%STATEBIN2FRAME Convert sleep score time bins to frame indices
%   Usage: idx = statebin2frame(bin1, bin2, ..., target_taxis)
%   Inputs:
%       binN: Nx2 matrix of [start_time, end_time]
%       target_taxis: vector of time points

if nargin < 2
    error('At least one timebin and target_taxis are required.');
end

target_taxis = varargin{end};

% Combine all timebins
sleepscore_timebin = [];
for k = 1:nargin-1
    sleepscore_timebin = [sleepscore_timebin; varargin{k}];
end

num_awake = size(sleepscore_timebin, 1);
time_idx = [];

for i = 1:num_awake
    [~, start_loc] = min(abs(target_taxis - sleepscore_timebin(i,1)));
    [~, end_loc] = min(abs(target_taxis - sleepscore_timebin(i,2)));
    time_idx = [time_idx, start_loc:end_loc];
end

% Ensure unique sorted indices
time_idx = unique(time_idx);
end

