function [pivtable, pivarray] = pivevent_assemble(piv_run, profile, dirrow, opt)
%PIVEVENT_ASSEMBLE  One session's PIV products as the two variables pivevent.mat holds.
%   Row of pivtable = one event. The cell blocks ride as multi-dimensional table
%   variables, so groupsummary averages a whole n_bin x n_wedge block at once and
%   the long form never has to be stored.
%
%   Once this has run the stack is not needed again for anything downstream of
%   the correlation. What still forces a re-read is written down in
%   design/PIVEVENT_DESIGN.md.
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
%      mean_ch1  H x W float, recording mean of channel 1. sd_ch1 likewise, and
%      sd_ch1    the same pair for channel 2. These four close the pixel-level
%      mean_ch2  questions that would otherwise reopen the stack
%      sd_ch2
%      pivparam  struct, every stage's settings. Stored, never read here
% OUT  pivtable  N x 28 table. pivparam sits in Properties.UserData, which is
%                what makes the bin columns readable -- a bin index means
%                nothing without bin_edges_um, and UserData cannot be separated
%                from the table the way a second variable can
%
%   ROW = (event, segment). Today every event is correlated as ONE displacement
%   from its start endpoint to its end endpoint, so seg_idx is 1 on every row and
%   seg_from / seg_to equal from / to. A field measured inside an event would be
%   several rows of the same event_idx with different seg_to, which needs no
%   change to this schema -- see CLAUDE_LOG.md for why the row unit is not the
%   event itself
%      pivarray  struct, the eight arrays. Saved as its own top-level variable
%                so load(f, 'pivtable') does not drag 15 MB with it
%
% UNIT  a suffix names a unit; no suffix is dimensionless or an index. dD_px is
%       the point of the rule -- the diameter change is pixels and has been
%       labelled microns twice, see CLAUDE_LOG.md
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

% Rows are matched by event id. piv_polar_events appends in piv_run's order, so
% position would work today, but the id is what the two records actually agree on
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

% 2. Per-row scalars. from / to are the event; seg_from / seg_to are the two
%    instants this row's displacement actually spans, which is the whole event
%    while an event is correlated as a single field
event_idx = zeros(n_event, 1);
seg_idx   = ones(n_event, 1);
state     = strings(n_event, 1);
polarity  = strings(n_event, 1);
bout      = zeros(n_event, 1);
gated     = false(n_event, 1);
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
    gated(k)     = row.gated;
    from(k)      = run.from;
    to(k)        = run.to;
    seg_from(k)  = run.from;
    seg_to(k)    = run.to;
    rise_s(k)    = run.rise_s;
    nfr(k)       = run.nfr;
    sgn(k)       = row.sgn;
    dD_px(k)     = row.dD;
end

% 3. Cell blocks. n_event x n_bin x n_wedge, the first dimension being the one
%    table requires to match its height. groupsummary then averages the block
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
pivtable = [keys, table(event_idx, seg_idx, state, polarity, bout, gated, ...
    from, to, seg_from, seg_to, rise_s, nfr, sgn, dD_px)];
for b = 1:numel(block_names)
    pivtable.(block_names(b)) = blocks.(block_names(b));
end
pivtable.first_bin = first_bin;
pivtable.reach_bin = reach_bin;

% 5. What the columns mean. The names carry the units for a reader; these slots
%    carry them for code, and the description says what the block's axes are
pivtable.Properties.UserData = opt.pivparam;
pivtable.Properties.VariableUnits = pivevent_units(pivtable.Properties.VariableNames);
pivtable.Properties.VariableDescriptions = ...
    pivevent_descriptions(pivtable.Properties.VariableNames, n_bin, n_wedge);

% 6. The arrays. xyuv is the record; everything in the table above was derived
%    from it, so a question the table cannot answer reopens here and not the stack
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
        'gated',     'whether uv came through the PIV gates'; ...
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
