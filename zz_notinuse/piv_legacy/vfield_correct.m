function uv_out = vfield_correct(uv_matrix, global_sub)
%VFIELD_CORRECT  Remove bulk motion from a PIV displacement field.
%   Subtracts one global median vector (whole-field drift, tissue motion, stage
%   translation) so only local/relative flow remains. The component index is the
%   LAST dimension: (:,:,...,1) = u, (:,:,...,2) = v. NaN off the sparse PIV grid
%   is allowed and stays NaN.
%
% IN   uv_matrix   H x W x 2 or H x W x N x 2 displacement field
%      global_sub  subtract the global median vector, default true
% OUT  uv_out      same size as uv_matrix, background subtracted

arguments
    uv_matrix  {mustBeNumeric}
    global_sub (1,1) logical = true
end

% 0. Pass the field through unchanged when subtraction is off
uv_out = uv_matrix;
if ~global_sub
    return;
end

% 1. Index the component dimension, so H x W x 2 and H x W x N x 2 both work
cdim = ndims(uv_matrix);
idx  = repmat({':'}, 1, cdim);

% 2. Split the two components
idx{cdim} = 1; u = uv_matrix(idx{:});
idx{cdim} = 2; v = uv_matrix(idx{:});

% 3. Background vector: median over all valid vectors
u_bg = median(u(:), 'omitnan');  if isnan(u_bg), u_bg = 0; end
v_bg = median(v(:), 'omitnan');  if isnan(v_bg), v_bg = 0; end

% 4. Write the corrected components back
idx{cdim} = 1; uv_out(idx{:}) = u - u_bg;
idx{cdim} = 2; uv_out(idx{:}) = v - v_bg;
end
