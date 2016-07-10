function [newBeginning, newEnd] = elongateDnaBackbone(PixelIdxList, bwThick)
    a = PixelIdxList;
    imgSize = size(bwThick);
    % for any entry, check in bwImgThickDna image where there 
    % are connected pixels around current pixel at idx with such a mask:
    % 1      1      1
    % 1     idx     1
    % 1      1      1
    % idx:  current pixel index
    % 1:    check in bwImgThickDna image whether this pixel is 1 or 0
    mask2 = [-imgSize(1)-1, -imgSize(1), -imgSize(1)+1, -1, ...
              1, imgSize(1)-1, imgSize(1), imgSize(1)+1 ];
    newBeginning = [];
    newEnd = [];
    % only if fragment has at least two pixels:
    if ~(size(a,1) < 2)
        % DNA fragment backbone now looks like this:
        % both end pixels only have one neighbor, all others only have 2
        % "endPixel1 - neighbor - ... - neighbor - endPixel2"
        endPixel1 = a(1);
        endPixel2 = a(size(a,1)); 
        %% first, we check whether there are any pixels BEFORE endPixel1
        mask_end1 = mask2;
        % set those entries from 8-neighborhood mask that are out of bounds
        % to zero
        mask_end1(endPixel1+mask_end1 < 1 | mask_end1+endPixel1 > imgSize(1)*imgSize(2) ) = 0;
        % now, we check the location of the center pixel's neighbor in the
        % 8-neighborhood to get the direction we have to check for "new"
        % pixels. It could look like this
        %       0      0      1
        %       0     idx     0
        %       0      0      0
        nextNeighbour1 = mask_end1+endPixel1 == ones(size(mask_end1))*a(2);
        % we invert the result:
        %       0      0      0
        %       0     idx     0
        %       1      0      0
        % ... and we therefore now know, how many pixel-indices we have to
        % subtract in order to get to the next pixel:
        %       0      0      idx
        %       0      1      0
        %       1      0      0
        % in this example, we would have to always substract imgSize(1)-1
        % many
        nextNeighbour1 = fliplr(nextNeighbour1);
        nextNeighbourDifference = mask_end1(find(nextNeighbour1));
        % if we end up at an out-of-bounds-index, the difference is 0 and
        % we don't have anything to do at this end
        if( nextNeighbourDifference ~= 0 )
            nextNeighbourIdx = endPixel1+nextNeighbourDifference;
            % we repeat the above steps as long as there are white pixels in
            % the direction we are going (in bwImgThickDna)
            while (bwThick(nextNeighbourIdx))
                newBeginning = [nextNeighbourIdx; newBeginning];
                oldIdx = nextNeighbourIdx;
                nextNeighbourIdx = nextNeighbourIdx+nextNeighbourDifference;
                if ( (nextNeighbourIdx < 1) ||                      ... % outside left side of image
                     (nextNeighbourIdx > imgSize(1)*imgSize(2)) ||  ... % outside right side of image
                     (mod(oldIdx, imgSize(1)) == 0 && mod(nextNeighbourIdx, imgSize(1)) == 1) ||  ... % outside lower side of image
                     (mod(oldIdx, imgSize(1)) == 1 && mod(nextNeighbourIdx, imgSize(1)) == 0) )       % outside upper side of image
                    break;
                end
            end
        end
        %% next, we check whether there are any pixels AFTER endPixel2
        % we do the same as above for endPixel1
        mask_end2 = mask2;
        mask_end2(endPixel2+mask_end2 < 1 | ...
            mask_end2+endPixel2 > imgSize(1)*imgSize(2) ) = 0;
        nextNeighbour2 = mask_end2+endPixel2 == ...
            ones(size(mask_end2))*a(size(a,1)-1);
        nextNeighbour2 = fliplr(nextNeighbour2);
        nextNeighbourDifference = mask_end2(find(nextNeighbour2));
        if( nextNeighbourDifference ~= 0)
            nextNeighbourIdx = endPixel2+nextNeighbourDifference;
            while (bwThick(nextNeighbourIdx))
                newEnd = [nextNeighbourIdx; newEnd];
                oldIdx = nextNeighbourIdx;
                nextNeighbourIdx = nextNeighbourIdx+nextNeighbourDifference;
                if ( (nextNeighbourIdx < 1) ||                      ... % outside left side of image
                     (nextNeighbourIdx > imgSize(1)*imgSize(2)) ||  ... % outside right side of image
                     (mod(oldIdx, imgSize(1)) == 0 && mod(nextNeighbourIdx, imgSize(1)) == 1) ||  ... % outside lower side of image
                     (mod(oldIdx, imgSize(1)) == 1 && mod(nextNeighbourIdx, imgSize(1)) == 0) )       % outside upper side of image
                    break;
                end
            end
        end
    end
end