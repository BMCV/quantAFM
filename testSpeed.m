function [i,gr, singlePath] = testSpeed(PixelIdxList, bwImg)
    len = size(PixelIdxList,2);
    res = [];
    for i= 1:len
        [gr, singlePath] = makeGraph(PixelIdxList{i}, bwImg);
    end