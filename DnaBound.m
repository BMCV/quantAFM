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
                 image, position, hasNucleus, nukleo)
            if nargin > 0
                dnaObj.connectedThick = connectedThick;
                dnaObj.bwImage = detail_thickDna;
                dnaObj.image = image;
                dnaObj.position = position;
                dnaObj.hasNucleus = hasNucleus;
%                 delete all other objects on the subImage besides the Dna
%                 strand
                dnaObj.isValid = 1;
                CC = bwconncomp(detail_thickDna);
                numPixels = cellfun(@numel, CC.PixelIdxList);
                [largestObj, idx] = max(numPixels);
                dnaObj.bwImage = false(size(dnaObj.bwImage));
                dnaObj.bwImage(CC.PixelIdxList{idx}) = true;
%%               Thin the subImage with Matlab function
                dnaObj.bwImageThinned = bwmorph(dnaObj.bwImage,'thin',Inf);
                
%               Calculate connected Components PixelIdxList
                connectedThinned = bwconncomp(dnaObj.bwImageThinned);
                if(~isempty(connectedThinned.PixelIdxList))
                    dnaObj.connectedThinned = connectedThinned.PixelIdxList{1};
                end
                dnaObj.sizeImg = size(detail_thickDna);
               
            end
%             set Reference to attached Nukeo
            if nargin > 5
                dnaObj.attachedNukleo = nukleo;
%                 nukleo.attachedTo = dnaObj;
            end
         end
    end
    
end

