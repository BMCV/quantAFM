
classdef nukleo
    %NUKLEO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        rad double;
        intensity double;   % max or mean
        volume double;
        position double = zeros(1,2);
        armRatio double; % Short / Long
        attachedTo DnaBound;  
    end
    
    methods
        function nukleoObj = nukleo(rad, intensity, volume, position, armRatio, attachedTo)
            if nargin > 0
                nukleoObj.rad = rad;
                nukleoObj.intensity = intensity;
                nukleoObj.volume = volume;
                nukleoObj.position = position;
                nukleoObj.armRatio = armRatio;
                
            end
            if nargon > 5
                nukleoObj.attachedTo = attachedTo;
                attachedTo.nukleo = nukleoObj;
            end
        end
        
       
    end
    
end

