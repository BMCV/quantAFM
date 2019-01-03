function [ adjMat ] = createAdjMat( picture, pixelIdxList )
% createAdjMat creates the AdjacencyMatrix for an object
% it needs the pixellist of the object and the picture on which it is
% defined with a line of 0s around the picture

% depending in which direction the next pixel is the entry in the matrix
% differs as followed
% 1 8 7
% 2 x 6
% 3 4 5

pixels = pixelIdxList;
[height, width] = size(picture);
adjMat = zeros(length(pixels), length(pixels));

for i = 1:length(pixels)
    % connected to the upper left pixel?
    if pixels(pixels == pixels(i) - height - 1)
        adjMat(find(pixels == pixels(i) - height - 1), i) = 1;
    end
    % connected to the left pixel?
    if pixels(pixels == pixels(i) - height)
        adjMat(find(pixels == pixels(i) - height), i) = 2;
    end
    % connected to the lower left pixel?
    if pixels(pixels == pixels(i) - height + 1)
        adjMat(find(pixels == pixels(i) - height + 1), i) = 3;
    end
    % connected to the upper pixel?
    if pixels(pixels == pixels(i) - 1)
        adjMat(find(pixels == pixels(i) - 1), i) = 8;
    end
    % connected to the lower pixel?
    if pixels(pixels == pixels(i) + 1)
        adjMat(find(pixels == pixels(i) + 1), i) = 4;
    end
    % connected to the upper right pixel?
    if pixels(pixels == pixels(i) + height - 1)
        adjMat(find(pixels == pixels(i) + height - 1), i) = 7;
    end
    % connected to the right pixel?
    if pixels(pixels == pixels(i) + height)
        adjMat(find(pixels == pixels(i) + height), i) = 6;
    end
    % connected to the lower right pixel?
    if pixels(pixels == pixels(i) + height + 1)
        adjMat(find(pixels == pixels(i) + height + 1), i) = 5;
    end
end

end

