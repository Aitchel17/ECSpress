function binary_array = statearray2binary(statecode,behavState)
%STATEARRAY2BINARY Summary of this function goes here
%   Detailed explanation goes here
    binary_array = double(behavState == statecode);
end

