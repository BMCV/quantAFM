classdef DnaBound < DNA
    %DNABOUND Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        attachedNukleo nukleo;
        angle double;
    end
    
    methods
         function dnaObj = DNA(connected, region, type, nukleo)
            if nargin > 0
                dnaObj.conntected = connected;
                dnaObj.region = region;
                dnaObj.position = region.Centroid;
                if max(strcmp(type,{normal, mutant}))==1
                   dnaObj.type = type;
                else
                   disp([type,' is not a valid DNA-type, please use: normal or mutant']);
                end
               
            end
            if nargin > 4
                dnaObj.nukleo = nukleo;
                nukleo.attachedTo = dnaObj;
            end
         end
    end
    
end

