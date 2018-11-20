function done = fusedImages( imgObj, imgName, os )
%FUSEDIMAGES creating additional overlay images of the output
%   Detailed explanation goes here
    thickfuse = imfuse(imgObj.rawImage, imgObj.bwImgThickDna);
    thinfuse = imfuse(imgObj.rawImage, imgObj.bwImgThinnedRemoved);
    bothfuse = imfuse(thickfuse, imgObj.bwImgThinnedRemoved);
    if (os == 0)
		imwrite(thickfuse, ['..\pictures\fused_images\' 'thickfuse_' imgName ]);
		imwrite(thinfuse, ['..\pictures\fused_images\' 'thinfuse_' imgName ]);
		imwrite(bothfuse, ['..\pictures\fused_images\' 'bothfuse_' imgName ]);    
    else
		imwrite(thickfuse, ['../pictures/fused_images/' 'thickfuse_' imgName ]);
		imwrite(thinfuse, ['../pictures/fused_images/' 'thinfuse_' imgName ]);
		imwrite(bothfuse, ['../pictures/fused_images/' 'bothfuse_' imgName ]); 
    end
    done = 'done';
end

