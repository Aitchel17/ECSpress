function xyuv_out = piv_blank(xyuv, keep)
%PIV_BLANK  NaN the displacement of the windows a mask excludes, keep the grid.
%   The sparse-to-sparse sibling of piv_stamp: same verdict, applied without
%   leaving the interrogation grid. piv_stamp writes the kept windows onto the
%   frame; this one blanks the rest where they sit.
%
% IN   xyuv      ny x nx x 4 float px  (:,:,1:2) = [x y] window centre,
%                                      (:,:,3:4) = [u v] displacement
%      keep      ny x nx logical       which windows survive
% OUT  xyuv_out  the same array with u and v NaN outside keep
%
% UNIT  whatever xyuv(:,:,3:4) was in. Nothing is scaled here
%   why  the COORDINATES stay. An array that dropped them could not tell a window
%        that was measured and rejected from one the grid never had, and both the
%        gate audit and the ungated comparison need that difference
%   why  one file rather than a local in each producer : analysis_pivensemble and
%        piv_run_events both apply a keep mask to a grid, and two copies of four
%        lines drift in exactly the way a gated and an ungated field must not

arguments
    xyuv (:,:,4) double
    keep (:,:) logical
end
if ~isequal(size(keep), [size(xyuv, 1), size(xyuv, 2)])
    error('piv_blank:keepShape', ...
        'keep is %s, the grid is %d x %d.', ...
        mat2str(size(keep)), size(xyuv, 1), size(xyuv, 2));
end

xyuv_out = xyuv;
u = xyuv_out(:,:,3);
v = xyuv_out(:,:,4);
u(~keep) = NaN;
v(~keep) = NaN;
xyuv_out(:,:,3) = u;
xyuv_out(:,:,4) = v;
end
