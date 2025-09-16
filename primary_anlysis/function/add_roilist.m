function roilist = add_roilist(loaded_data)
%ADD_ROILIST Summary of this function goes here
%   Detailed explanation goes here
    if isfield(loaded_data, 'roilist')
        roilist = loaded_data.roilist;
   
    else
        disp('roilist does not exist')
        roilist = struct();
    end
end

