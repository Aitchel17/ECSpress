function output_struct = secondary_afterproceesing(secondary_struct)
%SECONDARY_AFTERPROCEESING Summary of this function goes here
    % 1. Does total PVS decrease durin dilation?
    % 2. Calculate linear fit between PVS and BV relation
    % 3. Define dynamic and static PVS

%   Detailed explanation goes here
clean_ssidx_list = [];
for ssidx=1:size(secondary_struct,2)
    % does the total pvs decrease while vessel increase --> if not its not
    % perivascular space
    tmp.minpvschanges =  abs(min(secondary_struct(ssidx).heatdata(1).modespvs,[],'all'));
    tmp.maxpvschanges =  abs(max(secondary_struct(ssidx).heatdata(1).modespvs,[],'all'));
    if tmp.minpvschanges > tmp.maxpvschanges
        clean_ssidx_list = [clean_ssidx_list, ssidx];
        for idx = 1:3
           [b, ~] = robustfit(secondary_struct(ssidx).heatdata(idx).x_centers_aligned,...
                   secondary_struct(ssidx).heatdata(idx).modespvs);
           secondary_struct(ssidx).heatdata(idx).slope = b(2);
           secondary_struct(ssidx).heatdata(idx).intercept = b(1);
           midpoint = floor(size(secondary_struct(ssidx).heatdata(idx).x_centers_aligned,2)/2);

           % first half
           try
               bvdiam_change  = secondary_struct(ssidx).heatdata(idx).x_centers_aligned(1:midpoint);
               pvsthickness_change = secondary_struct(ssidx).heatdata(idx).modespvs(1:midpoint);
               [b1, ~] = robustfit(bvdiam_change,pvsthickness_change);
               secondary_struct(ssidx).heatdata(idx).slope_firsthalf = b1(2);
               secondary_struct(ssidx).heatdata(idx).intercept_firsthalf = b1(1);
               disp(midpoint)

  
           catch
               % disp([ssidx idx])
           end
            try
               bvdiam_change  = secondary_struct(ssidx).heatdata(idx).x_centers_aligned(midpoint:end);
               pvsthickness_change = secondary_struct(ssidx).heatdata(idx).modespvs(midpoint:end);
               [b2, ~] = robustfit(bvdiam_change,pvsthickness_change);
               secondary_struct(ssidx).heatdata(idx).slope_secondhalf = b2(2);
               secondary_struct(ssidx).heatdata(idx).intercept_secondhalf = b2(1);
                disp(midpoint)
            catch
            end
        end

            % swap row 1 and row 2, so row 1 to be always dynamic pvs
        if secondary_struct(ssidx).heatdata(1).slope > secondary_struct(ssidx).heatdata(2).slope
            disp(secondary_struct(ssidx).heatdata(1).slope)
            secondary_struct(ssidx).heatdata = secondary_struct(ssidx).heatdata([2,1,3]);
        end
    else
        % disp(ssidx)
    end


end
output_struct = secondary_struct(clean_ssidx_list);

end

