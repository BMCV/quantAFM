function [i,gr, singlePath] = testSpeed(PixelIdxList, bwImg, bwThick)
    len = size(PixelIdxList,2);
    res = zeros(size(bwImg));
    for i= 1:len
        [gr, singlePath] = makeGraph(PixelIdxList{i}, bwImg, bwThick);
        res(singlePath) = 1;
    end
    imshow(imfuse(bwThick, res));
    imwrite(imfuse(bwThick, res) , '../pictures/DNA_spine_thickDna_overlay.tif');