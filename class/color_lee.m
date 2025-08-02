classdef color_lee
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        red_g = []
        green_g = []
        blue_g = []
        magenta_g = []
        cyan_g = []
        yellow_g = []
        gray_g = []
    end
    
    methods
        function obj = color_lee()
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.red_g =    [linspace(0,1,256)', zeros(256,1), zeros(256,1)];
            obj.green_g =    [zeros(256,1), linspace(0,1,256)', zeros(256,1)];
            obj.blue_g =    [zeros(256,1), zeros(256,1), linspace(0,1,256)'];
            obj.magenta_g =    [linspace(0,1,256)', zeros(256,1), linspace(0,1,256)'];
            obj.cyan_g =    [zeros(256,1),linspace(0,1,256)', linspace(0,1,256)'];
            obj.yellow_g =    [linspace(0,1,256)', linspace(0,1,256)',zeros(256,1)];
            obj.gray_g =    [linspace(0,1,256)', linspace(0,1,256)', linspace(0,1,256)'];
        end
    end
end

