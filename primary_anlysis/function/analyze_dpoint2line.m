function [orthogonal_vertice, dist] = analyze_dpoint2line(point_vertice,line_vertices)
%ANALYZE_DPOINT2LINE Summary of this function goes here
%   Detailed explanation goes here
    lstart = line_vertices(1,:);
    lend = line_vertices(2,:);
    origin = point_vertice;
    start2end = lend-lstart;
    start2point = origin - lstart;
    projected_length = dot(start2end,start2point);
    totallength = dot(start2end,start2end);  
    ratio = projected_length/totallength;
    orthogonal_vertice =  lstart + ratio*start2end;
end

