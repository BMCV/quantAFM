classdef DnaFree < DNA
    %DNAFREE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function dnaObj = DnaFree(connectedThick,subImage, position)
            if nargin > 0
                dnaObj.position = position;
                dnaObj.connectedThick = connectedThick;
                dnaObj.bwImage = subImage;
         
                dnaObj.type = 'free';
                                
%                 delete all other objects on the subImage besides the Dna
%                 strand
                CC = bwconncomp(subImage);
                numPixels = cellfun(@numel, CC.PixelIdxList);
                [largestObj, idx] = max(numPixels);
                dnaObj.bwImage = false(size(dnaObj.bwImage));
                dnaObj.bwImage(CC.PixelIdxList{idx}) = true;
%               Thinn the subImage
                dnaObj.bwImageThinned = bwmorph(dnaObj.bwImage,'thin',Inf);
%               Calculate connected Components PixelIdxList
                connectedThinned = bwconncomp(dnaObj.bwImageThinned);
                dnaObj.connectedThinned = connectedThinned.PixelIdxList{1};
                dnaObj.sizeImg = size(subImage);
               
            end
        end
    end
    
end

