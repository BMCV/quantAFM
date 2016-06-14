classdef DnaFree < DNA
    %DNAFREE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function dnaObj = DnaFree(connectedThick,connectedThinned, position)
            if nargin > 0
                dnaObj.position = position;
                dnaObj.connectedThick = connectedThick;
                dnaObj.connectedThinned = connectedThinned;
         
                dnaObj.type = 'free';
             
                 %   disp([type,' is not a valid DNA-type, please use: free, normal or mutant']);
               
            end
        end
    end
    
end

