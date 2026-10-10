function time_table = statebin2timetable(varargin)
%STATEBIN2BOUTSEC Summary of this function goes here
% Input: time bins dim: [startbin, endbin]
% Output: combined and segmented time table
%   Detailed explanation goes here

    combined_bintable = cat(1,varargin{:});
    %%
    sorted_bintable = sort(combined_bintable,1);

    nbin = length(sorted_bintable);
    nbin_col = 1;
    if ~isempty(sorted_bintable)
        time_table = sorted_bintable(1,:);
    
        for nidx = 2:nbin
            if time_table(nbin_col,2) == sorted_bintable(nidx,1)
                time_table(nbin_col,2) = sorted_bintable(nidx,2);
            else
                nbin_col = nbin_col+1;
                time_table(nbin_col,:) = sorted_bintable(nidx,:);
            end
        end
    else
        time_table = [];
    end
end

