
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
%         Number of DNA strand it is attached to
        number;
    end
    
    methods
        function nukleoObj = nukleo(center, rad, number)
            if nargin > 0
                nukleoObj.rad = rad;
                nukleoObj.center = center;
                nukleoObj.number = number;
                
            end
            
        end
        
       
    end
    
end

