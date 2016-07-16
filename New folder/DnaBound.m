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
                 position, hasNucleus, nukleo)
            if nargin > 0
                dnaObj.connectedThick = connectedThick;
                dnaObj.bwImage = detail_thickDna;

                dnaObj.position = position;
                dnaObj.hasNucleus = hasNucleus;
%                 delete all other objects on the subImage besides the Dna
%                 strand
                CC = bwconncomp(detail_thickDna);
                numPixels = cellfun(@numel, CC.PixelIdxList);
                [largestObj, idx] = max(numPixels);
                dnaObj.bwImage = false(size(dnaObj.bwImage));
                dnaObj.bwImage(CC.PixelIdxList{idx}) = true;
%%               Thin the subImage with Matlab function
%                dnaObj.bwImageThinned = bwmorph(dnaObj.bwImage,'thin',Inf);
%%               ... using user-written implementation
                dnaObj.bwImageThinned = thinningZhangSuen(dnaObj.bwImage);
                % visualize
                % dna1 = bwmorph(dnaObj.bwImage,'thin',Inf);dna2 = thinningZhangSuen(dnaObj.bwImage);figure; subplot(3,1,1); imshow(dna1); subplot(3,1,2); imshow(dna2); subplot(3,1,3); imshow(dnaObj.bwImage);

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

