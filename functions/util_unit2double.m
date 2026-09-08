function value = util_unit2double(text)
%UTIL_UNIT2DOUBLE  The leading number of a unit-suffixed value, as a double.
%   Caller: centralize_driftstate
%   The same function as 00_mdfExtractor/mdfExtractor/function/util_unit2double,
%   taking a column as well as a scalar so a table variable goes through in one call.
%
% IN   text   char | string | cell | double   '70um', '0.285 µm', a number, or an
%                                             N x 1 column of those
% OUT  value  N x 1 double                    NaN where no number leads the text
    if isa(text, 'double')
        value = text;
        return
    end
    text = string(text);
    value = nan(size(text));
    for k = 1:numel(text)
        head = sscanf(text(k), '%f', 1);
        if ~isempty(head)
            value(k) = head;
        end
    end
end
