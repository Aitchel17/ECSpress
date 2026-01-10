classdef bleach < roi
    %BLEACH Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function obj = bleach(reference_stack,roimod,mdfExtractLoader_instance)
            arguments
                reference_stack (:,:,:) {mustBeNumeric}
                roimod (1,:) char {mustBeMember(roimod, ["polygon", "rectangle"])}
                mdfExtractLoader_instance mdfExtractLoader = []
            end

           obj@roi(reference_stack,roimod,mdfExtractLoader_instance);
        end
        
        function data = decay(obj,stackfieldname)
            data = obj.data;
            % field name construction
            fname_ave = ['ave_', stackfieldname];
            fname_norm = ['norm_', stackfieldname];
            fname_exponential1 = ['exponential1_', stackfieldname];
            fname_exponential2 = ['exponential2_', stackfieldname];

            data.(fname_ave) = squeeze(mean(obj.stacks.(stackfieldname),[1,2],'omitmissing'));
            data.(fname_norm) = data.(fname_ave)/data.(fname_ave)(1);
            data.(fname_exponential1) = log(data.(fname_ave));
            data.(fname_exponential1) = data.(fname_exponential1)/data.(fname_exponential1)(1);
            data.(fname_exponential2) = log(log(data.(fname_ave)));
            data.(fname_exponential2) = data.(fname_exponential2)/data.(fname_exponential2)(1);
        end
    end
end

