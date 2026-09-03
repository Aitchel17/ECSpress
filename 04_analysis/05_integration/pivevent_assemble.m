function [pivtable, pivarray] = pivevent_assemble(piv_run, profile, dirrow, opt)
%PIVEVENT_ASSEMBLE  One session's PIV products as the two variables pivevent.mat holds.
%   Row of pivtable = one event; the cell blocks ride as multi-dimensional table
%   variables, so groupsummary averages a whole n_bin x n_wedge block at once.
%   What still forces a re-read of the stack : design/PIVEVENT_DESIGN.md
%
% IN   piv_run   1 x N struct, piv_run_events' output. Event metadata and the
%                two displacement grids
%      profile   vfield_profile holding one row per event. Rows are matched to
%                piv_run BY event_idx, not by position
%      dirrow    1 x 1 table, the session's row of the directory table. Its key
%                columns are copied verbatim so a join against the FWHM table
%                needs no translation
%      coremask  H x W bool, the vessel at maximum dilation
%      exclmask  H x W bool, true = pixel excluded from PIV
%      mean_ch1  H x W float, recording mean of channel 1; sd_ch1, mean_ch2,
%      sd_ch1    sd_ch2 likewise
%      mean_ch2
%      sd_ch2
%      pivparam  struct, every stage's settings. Stored, never read here
% OUT  pivtable  N x 28 table. pivparam sits in Properties.UserData, so a bin
%                index can be read against bin_edges_um without a second variable
%      pivarray  struct, the eight arrays, saved as its own top-level variable
%
%   ROW = (event, segment). Today seg_idx is 1 on every row and seg_from / seg_to
%   equal from / to; a field measured inside an event would add rows of the same
%   event_idx. Why the row unit is not the event : see CLAUDE_LOG.md
% UNIT  a suffix names a unit; no suffix = dimensionless or an index. dD_px has
%       been labelled microns twice, see CLAUDE_LOG.md
%
%   see also PIV_RUN_EVENTS, PIV_POLAR_EVENTS, VFIELD_PROFILE

arguments
    piv_run      (1,:) struct
    profile      (1,1) vfield_profile
    dirrow       (1,1) table
    opt.coremask logical
    opt.exclmask logical
    opt.mean_ch1 double
    opt.sd_ch1   double
    opt.mean_ch2 double
    opt.sd_ch2   double
    opt.pivparam (1,1) struct
end

n_event = numel(piv_run);
n_bin   = profile.n_bin;
n_wedge = profile.n_wedge;

% rows are matched by event id, not by position
row_id  = [profile.rows.event_idx];
run_id  = [piv_run.id];
[found, row_of] = ismember(run_id, row_id);
if ~all(found)
    error('pivevent_assemble:noProfileRow', ...
        'piv_run has event %d, which profile has no row for.', run_id(find(~found, 1)));
end

% 1. Key columns, copied from the directory row without translation
key_names = ["MouseID", "Date", "SessionType", "SessionID", "Directory", "VesselID"];
keys      = table();
for name = key_names
    value = dirrow.(name);
    keys.(name) = repmat(string(value), n_event, 1);
end

% 2. Per-row scalars. seg_from / seg_to are the instants this row spans : the whole event today
event_idx = zeros(n_event, 1);
seg_idx   = ones(n_event, 1);
state     = strings(n_event, 1);
polarity  = strings(n_event, 1);
bout      = zeros(n_event, 1);
gate_name = strings(n_event, 1);
from      = zeros(n_event, 1);
to        = zeros(n_event, 1);
seg_from  = zeros(n_event, 1);
seg_to    = zeros(n_event, 1);
rise_s    = zeros(n_event, 1);
nfr       = zeros(n_event, 1);
sgn       = zeros(n_event, 1);
dD_px     = zeros(n_event, 1);

for k = 1:n_event
    run = piv_run(k);
    row = profile.rows(row_of(k));
    event_idx(k) = run.id;
    state(k)     = string(run.state);
    polarity(k)  = string(run.pol);
    bout(k)      = run.bout;
    gate_name(k) = row.gate_name;
    from(k)      = run.from;
    to(k)        = run.to;
    seg_from(k)  = run.from;
    seg_to(k)    = run.to;
    rise_s(k)    = run.rise_s;
    nfr(k)       = run.nfr;
    sgn(k)       = row.sgn;
    dD_px(k)     = row.dD;
end

% 3. Cell blocks, n_event x n_bin x n_wedge : the first dimension is the table height
block_names = ["divergence", "volume_out_um2", "disp_radial_um", "disp_tangential_um", ...
               "strain_radial", "strain_hoop", "strain_shear", "rotation", ...
               "n_tri", "area_um2"];
cell_names  = ["divergence", "volume_out",     "disp_radial",    "disp_tangential", ...
               "strain_radial", "strain_hoop", "strain_shear", "rotation", ...
               "n_tri", "area"];
blocks = struct();
for b = 1:numel(block_names)
    blocks.(block_names(b)) = nan(n_event, n_bin, n_wedge);
end
first_bin = nan(n_event, n_wedge);
reach_bin = nan(n_event, n_wedge);

for k = 1:n_event
    cells = profile.rows(row_of(k)).cells;
    for b = 1:numel(block_names)
        blocks.(block_names(b))(k, :, :) = cells.(cell_names(b));
    end
    first_bin(k, :) = cells.first_bin;
    reach_bin(k, :) = cells.reach_bin;
end

% 4. The table
pivtable = [keys, table(event_idx, seg_idx, state, polarity, bout, gate_name, ...
    from, to, seg_from, seg_to, rise_s, nfr, sgn, dD_px)];
for b = 1:numel(block_names)
    pivtable.(block_names(b)) = blocks.(block_names(b));
end
pivtable.first_bin = first_bin;
pivtable.reach_bin = reach_bin;

% 5. What the columns mean, for code : units and the block axes
pivtable.Properties.UserData = opt.pivparam;
pivtable.Properties.VariableUnits = pivevent_units(pivtable.Properties.VariableNames);
pivtable.Properties.VariableDescriptions = ...
    pivevent_descriptions(pivtable.Properties.VariableNames, n_bin, n_wedge);

% 6. The arrays. xyuv is the record everything above was derived from
pivarray = struct( ...
    'xyuv',         {{piv_run.xyuv}}, ...          % 1 x N cell, each ny x nx x 4 px
    'xyuv_ungated', {{piv_run.xyuv_ungated}}, ...  % 1 x N cell, before the gates
    'coremask',     opt.coremask, ...              % H x W bool
    'exclmask',     opt.exclmask, ...              % H x W bool
    'mean_ch1',     opt.mean_ch1, ...              % H x W float
    'sd_ch1',       opt.sd_ch1, ...                % H x W float
    'mean_ch2',     opt.mean_ch2, ...              % H x W float
    'sd_ch2',       opt.sd_ch2);                   % H x W float
end

% ---------------------------------------------------------------------------
function units = pivevent_units(names)
%PIVEVENT_UNITS  The unit slot for each column, empty where the value has none.
    known = { ...
        'rise_s',             's'; ...
        'dD_px',              'px'; ...
        'from',               'frame'; ...
        'to',                 'frame'; ...
        'seg_from',           'frame'; ...
        'seg_to',             'frame'; ...
        'volume_out_um2',     'um^2'; ...
        'disp_radial_um',     'um'; ...
        'disp_tangential_um', 'um'; ...
        'area_um2',           'um^2'};
    units = repmat({''}, 1, numel(names));
    for k = 1:numel(names)
        hit = strcmp(known(:,1), names{k});
        if any(hit)
            units{k} = known{hit, 2};
        end
    end
end

% ---------------------------------------------------------------------------
function desc = pivevent_descriptions(names, n_bin, n_wedge)
%PIVEVENT_DESCRIPTIONS  One line per column, saying what a block's axes are.
    block_note  = sprintf('height x %d bin x %d wedge, area-weighted per cell', ...
        n_bin, n_wedge);
    wedge_note  = sprintf('height x %d wedge, bin index', n_wedge);
    block_names = {'divergence', 'volume_out_um2', 'disp_radial_um', ...
        'disp_tangential_um', 'strain_radial', 'strain_hoop', 'strain_shear', ...
        'rotation', 'n_tri', 'area_um2'};
    known = { ...
        'event_idx', 'index into the eventlist analysis_event produced'; ...
        'seg_idx',   'which segment of the event this row measures. 1 while an event is one field'; ...
        'seg_from',  'the earlier instant this row''s displacement spans'; ...
        'seg_to',    'the later one. Equals to while seg_idx is always 1'; ...
        'state',     'transition label, or the state name for a quiet control'; ...
        'polarity',  'dilation / constriction / none. none = a control row'; ...
        'gate_name', 'the PIV gates uv came through, joined by +'; ...
        'nfr',       'frame pairs the ensemble averaged. NOT a scale factor'; ...
        'sgn',       'sign of dD_px, +1 or -1'; ...
        'dD_px',     'diameter change on the smoothed trace, PIXELS'; ...
        'first_bin', wedge_note; ...
        'reach_bin', wedge_note};
    desc = repmat({''}, 1, numel(names));
    for k = 1:numel(names)
        if any(strcmp(block_names, names{k}))
            desc{k} = block_note;
            continue
        end
        hit = strcmp(known(:,1), names{k});
        if any(hit)
            desc{k} = known{hit, 2};
        end
    end
end
