function hPatches = plot_sleep_patches(ax, state_data, varargin)
%PLOT_SLEEP_PATCHES Draw colored patches for sleep states on an axis
%   hPatches = plot_sleep_patches(ax, state_data)
%   hPatches = plot_sleep_patches(ax, state_data, t_axis)
%
%   Inputs:
%       ax: Target axes handle.
%       state_data: Struct where fields are state names (e.g., 'nrem', 'rem', 'awake')
%                   and values are Nx2 matrices of [Start, End].
%                   - If 't_axis' is provided, these are INDICES into t_axis.
%                   - If 't_axis' is NOT provided, these are TIME (seconds).
%       t_axis: (Optional) Time vector corresponding to indices in state_data.
%
%   Optional Name-Value pairs:
%       'YLim', [min max]   : Specific Y-limits for patches (default: ax.YLim).
%       'FaceAlpha', val    : Transparency (default 0.4).
%       'Tag', str          : Tag for patch objects (default 'SleepPatch').
%
%   Returns:
%       hPatches: Array of patch handles created.

% Parse Input
p = inputParser;
addRequired(p, 'ax', @(x) isgraphics(x, 'axes'));
addRequired(p, 'state_data', @isstruct);
addOptional(p, 't_axis', [], @isnumeric);
addParameter(p, 'YLim', [], @isnumeric);
addParameter(p, 'FaceAlpha', 0.4, @isnumeric);
addParameter(p, 'Tag', 'SleepPatch', @ischar);

parse(p, ax, state_data, varargin{:});

t_axis = p.Results.t_axis;
force_ylim = p.Results.YLim;
alpha_val = p.Results.FaceAlpha;
tag_val = p.Results.Tag;

% Define Colors (Internal Standard)
cols.nrem = [0.00 0.45 0.90];     % Blue
cols.rem = [0.85 0.10 0.10];      % Red
cols.awake = [0.10 0.10 0.10];    % Black
cols.drowsy = [0.20 0.65 0.20];   % Green
cols.uarousal = [0.93 0.49 0.19]; % Orange
cols.trans = [0.50 0.10 0.70];    % Purple

% Determine Y-Limits
if isempty(force_ylim)
    force_ylim = get(ax, 'YLim');
end
yLow = force_ylim(1);
yHigh = force_ylim(2);

% Hold axis
wasHeld = ishold(ax);
hold(ax, 'on');

state_names = fieldnames(state_data);
hPatches = gobjects(0);

for i = 1:length(state_names)
    sName = state_names{i};
    intervals = state_data.(sName);

    if isempty(intervals) || ~isnumeric(intervals) || size(intervals, 2) ~= 2
        continue;
    end

    % Determine Color based on name
    c = [0.5 0.5 0.5]; % Default Gray
    lowerName = lower(sName);

    if contains(lowerName, 'nrem')
        c = cols.nrem;
    elseif contains(lowerName, 'rem')
        c = cols.rem;
    elseif contains(lowerName, 'awake')
        c = cols.awake;
    elseif contains(lowerName, 'drowsy')
        c = cols.drowsy;
    elseif contains(lowerName, 'arousal')
        c = cols.uarousal;
    elseif contains(lowerName, 'trans')
        c = cols.trans;
    end

    % Iterate through bouts
    nBouts = size(intervals, 1);

    % Pre-calculate coordinates to avoid loop overhead for graphics?
    % Patch can take matrices, but usually separate patches are easier to manage/hit-test.
    % For performance with many patches, we can create one multipart patch or many patches.
    % Individual patches are better for transparency overlap (not an issue here usually)
    % and individual identification.

    for b = 1:nBouts
        sIdx = intervals(b, 1);
        eIdx = intervals(b, 2);

        % Convert to Time
        if ~isempty(t_axis)
            % Boundary Checks
            sIdx = max(1, min(sIdx, length(t_axis)));
            eIdx = max(1, min(eIdx, length(t_axis)));

            t1 = t_axis(sIdx);
            t2 = t_axis(eIdx);
        else
            t1 = sIdx;
            t2 = eIdx;
        end

        % Draw Patch
        hp = patch(ax, [t1 t2 t2 t1], [yLow yLow yHigh yHigh], c, ...
            'FaceAlpha', alpha_val, 'EdgeColor', 'none', ...
            'Tag', [tag_val '_' sName], 'HitTest', 'off');

        hPatches(end+1) = hp;
    end
end

if ~wasHeld
    hold(ax, 'off');
end

% Restore YLim (sometimes patch auto-scales)
ylim(ax, force_ylim);

end
