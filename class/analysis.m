classdef analysis
    %ANALYSIS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Property1
    end
    
    methods
        function obj = analysis(inputArg1,inputArg2)
            %ANALYSIS Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end

    properties (Constant)
            clist = struct(...
            'white', [1,1,1],...
            'black', [0,0,0],...
            'cyan', [0,1,1],...
            'magenta', [1,0,1],...
            'yellow', [1,1,0],...
            'red', [1,0,0],...
            'green', [0,1,0],...
            'blue', [0,0,1],...
            'orange', [1,0.5,0],...
            'lightgreen', [0.4,0.8,0.4],...
            'darkgreen', [0,0.5,0]);
    end
end

