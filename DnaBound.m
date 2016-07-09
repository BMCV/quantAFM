classdef DnaBound < DNA
    %DNABOUND Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
%          List of attached Nukleos Objects
        attachedNukleo;
        angle1 double;
        angle2 double;
    end
    
%     If attached Nukleo are provided(More than 4 inout arguments), Set them.
    methods
         function dnaObj = DnaBound(connectedThick,detail_thickDna, ...
                 position, type, nukleo)
            if nargin > 0
                dnaObj.connectedThick = connectedThick;
                dnaObj.bwImage = detail_thickDna;

                dnaObj.position = position;
                if max(strcmp(type,{'normal', 'mutant'}))==1
                   dnaObj.type = type;
                else
                   disp([type,' is not a valid DNA-type, please use: normal or mutant']);
                end
%                 delete all other objects on the subImage besides the Dna
%                 strand
                CC = bwconncomp(detail_thickDna);
                numPixels = cellfun(@numel, CC.PixelIdxList);
                [largestObj, idx] = max(numPixels);
                dnaObj.bwImage = false(size(dnaObj.bwImage));
                dnaObj.bwImage(CC.PixelIdxList{idx}) = true;
%%               Thin the subImage with Matlab function
                dnaObj.bwImageThinned = bwmorph(dnaObj.bwImage,'thin',Inf);
%%               ... using user-written implementation
%                dnaObj.bwImageThinned = thinningZhangSuen(dnaObj.bwImage);
                
%               Calculate connected Components PixelIdxList
                connectedThinned = bwconncomp(dnaObj.bwImageThinned);
                dnaObj.connectedThinned = connectedThinned.PixelIdxList{1};
                dnaObj.sizeImg = size(detail_thickDna);
               
            end
%             set Reference to attached Nukeo
            if nargin > 4
                dnaObj.attachedNukleo = nukleo;
%                 nukleo.attachedTo = dnaObj;
            end
         end
    end
    
end

