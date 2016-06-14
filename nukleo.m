
classdef nukleo
    %NUKLEO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
%         Radius of the nukleo fittedcircle
        rad double;
%         center of the fitted circle
        center = zeros(1,2);
        armRatio double; % Short / Long
        attachedTo;  
    end
    
    methods
        function nukleoObj = nukleo(rad, center, attachedTo)
            if nargin > 0
                nukleoObj.rad = rad;
                nukleoObj.center = center;
                
            end
            if nargin > 2
%                 nukleoObj.attachedTo = attachedTo;
%                 attachedTo.nukleo = nukleoObj;
            end
        end
        
       
    end
    
end

