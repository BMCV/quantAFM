classdef (Abstract) DNA
    %DNA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        number int;
        length double;
        
        position double = zeros(1,2) ; %center of mass
        type char;
        connected;
        region;
        
    end
    
    methods
    end
    
end


