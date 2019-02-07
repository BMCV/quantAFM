function done = fusedImages( imgObj, imgName, filepath)
%FUSEDIMAGES creating additional overlay images of the output
%   Detailed explanation goes here
    thickfuse = imfuse(imgObj.rawImage, imgObj.bwImgThickDna);
    thinfuse = imfuse(imgObj.rawImage, imgObj.bwImgThinnedRemoved);
    bothfuse = imfuse(thickfuse, imgObj.bwImgThinnedRemoved);
    imwrite(thickfuse, fullfile(filepath, 'fused_images', ['thickfuse_' imgName ]));
    imwrite(thinfuse, fullfile(filepath, 'fused_images', ['thinfuse_' imgName ]));
    imwrite(bothfuse, fullfile(filepath, 'fused_images', ['bothfuse_' imgName ])); 
    
    done = 'done';
end

