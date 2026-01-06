function [outputArg1,outputArg2] = sleepscore(sleepscore_inputdata)
%SLEEPSCORE Summary of this function goes here
% interface for determine the start and end time of sleep
% ECoG powerspectrum, Pupil diameter, EMG, and whisker movement will be
% ploted
% using slider bar, positioning the start and end bar, 

%   Detailed explanation goes here
% input data struct should be struct contains ECoG powerspectrum

Window = uifigure('Name', 'sleep scoring window');
MainLayout = uigridlayout(Window, [2,1]);  % Three panels (plot, uislider, buttons for sleep state)
MainLayout.RowHeight = {'1x', 100, 100};
MainLayout.ColumnWidth = {'1x'};

% images

end

