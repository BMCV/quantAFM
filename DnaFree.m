classdef DnaFree < DNA
    %DNAFREE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function dnaObj = DNA(connected, region)
            if nargin > 0
                dnaObj.position = region.Centroid;
                dnaObj.connected = connected;
                dnaObj.region = region;
                 
                dnaObj.type = 'free';
             
                 %   disp([type,' is not a valid DNA-type, please use: free, normal or mutant']);
               
            end
        end
    end
    
end

