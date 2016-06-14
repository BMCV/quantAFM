classdef DnaBound < DNA
    %DNABOUND Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
%          List of attached Nukleos Objects
        attachedNukleo;
        angle double;
    end
    
%     If attached Nukleo are provided(More than 4 inout arguments), Set them.
    methods
         function dnaObj = DnaBound(connectedThick,connectedThinned,position, type, nukleo)
            if nargin > 0
                dnaObj.connectedThick = connectedThick;
                dnaObj.connectedThinned = connectedThinned;

                dnaObj.position = position;
                if max(strcmp(type,{'normal', 'mutant'}))==1
                   dnaObj.type = type;
                else
                   disp([type,' is not a valid DNA-type, please use: normal or mutant']);
                end
               
            end
%             set Reference to attached Nukeo
            if nargin > 4
                dnaObj.attachedNukleo = nukleo;
%                 nukleo.attachedTo = dnaObj;
            end
         end
    end
    
end

