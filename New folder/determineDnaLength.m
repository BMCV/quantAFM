function [lengths, backboneIdxList] = determineDnaLength(PixelIdxList, bwImgThin, bwImgThick)
    len = size(PixelIdxList,2);
    res = zeros(size(bwImgThin)); %can be deleted later, only for visualization
    lengths = zeros(size(PixelIdxList));
    backboneIdxList = {};
    
%    imshow(res);
    %imshow(bwImgThin);
%    hold on;
    % for each DNA fragment in bwImgThin:
    for i= 1:len
        currPxlList = PixelIdxList{i};
        % if fragment too small, don't compute the rest
        if( size(currPxlList,1) < 3)
            continue
        end
        % create graph from its PixelIdxList and, from that, get the
        % fragment backbone
        [gr, singlePath] = getDnaBackbone(currPxlList, bwImgThin);
        % both end parts of the fragment were lost during the
        % thinning step, so get them back from the thickDNA fragment
        [newBeginning, newEnd] = ... 
            elongateDnaBackbone(singlePath, currPxlList, bwImgThick);
        % ... and add the respective pixels to the DNA backbone
        singlePath = [newBeginning; currPxlList(singlePath); newEnd];
        % calculate cubic splines for current backbone
        [row, col]= ind2sub(size(res),singlePath');
%        spline = cscvn([col; row]);
%        fnplt(spline, 2);
%        plot(col, row);
%        s = fnplt(spline);
        % calculate length of spline
        %[~, s1] = unique(round(s/1e-15),'rows','stable');
        %s2 = s(s1);
        %sqrt(sum(diff(s,[],1).^2,2))
        % the above is equal to: 
        %lengths(1,i)= pdist(s);
        
        %% calculate length with Kulpa Estimator
        % Length = 0.948*Ne + 1.343*No
        % with:     Ne - number of even pixels
        %           No - number of odd pixels
        % Here, even pixels are those that are, in their 8-neighourhood
        % connected along the even axes, and odd pixels are those connected
        % along the uneven axes.
        %   3   2   1        0   1   1      would have:
        %   4  pxl  0   =>   0   1   0  =>  2 odd pixels
        %   5   6   7        1   0   0      2 even pixel
        [row, col]= ind2sub(size(res),singlePath');
        evenOdds = zeros(size(col));
        % calculate differences between neighbouring entries in column and
        % in row indices, respectively. Even pixels should have difference
        % "0" between their row and col indices, resp.
        diff_col = diff(col');
        diff_row = diff(row');
        % now, create vector that has "1" for even and "0" for odd pixels
        evenOdds(find(~(diff_col))) = 1; 
        evenOdds(find(~(diff_row))) = 1;
        numberOfEvenPixels = sum(evenOdds);
        numberOfOddPixels = size(col, 2) - numberOfEvenPixels;
        % apply Kulpa Estimator
        fragmentLength = 0.948*numberOfEvenPixels + 1.343 * numberOfOddPixels;
        lengths(1,i) = fragmentLength;
        backboneIdxList{i} = singlePath;
        
        res(singlePath) = 1; %can be deleted later, only for visualization
    end
    'done'
%    imshow(imfuse(bwImgThick, res));
%    imwrite(imfuse(bwImgThick, res) , '../pictures/DNA_spine_thickDna_overlay.tif');
end

function [gr, singlePath ]= getDnaBackbone(pxlIdxList, bwImg)
    % a is shorter than pxlIdxList... :-)
    a = pxlIdxList;    
    imgSize = size(bwImg);
%    gr = diag(ones(1,size(a,1)));
    % future upper adjacency matrix
    gr = zeros(size(a,1),size(a,1));
    % for any entry, check in bwImgThinnedDna image where there 
    % are connected pixels around current pixel at idx with such a mask:
    % x      x      1
    % x     idx     1
    % x      1      1
    % x:    don't care, 
    % idx:  current pixel index
    % 1:    check in bwImgThinnedDna image whether this pixel is 1 or 0
    mask = [1, imgSize(1)-1, imgSize(1), imgSize(1)+1 ];
    % iterate over each entry in pixelIdxList
    for i = 1:size(a,1)
        tmp_mask = mask;
        % we have to check that we don't index out of bounds; so
        % respective mask values are removed
        tmp_mask(tmp_mask+a(i) > imgSize(1)*imgSize(2) ) = [];
        tmp_mask = tmp_mask + a(i);
        % apply mask and get indices where mask is true
        res = bwImg(tmp_mask).*(tmp_mask);
        % we need to exclude border cases; so, we test whether 
        % the found indices really occur in pixelIdxList 
        % ismemeber returns 0 and 1, respectively. So, we can
        % use its output as entry in adjacency matrix
        gr(i,:) = ismember(a, res(find(res)));
    end
    % alternative for for-loop ?:
    %c = bsxfun(@plus, a, mask);
    %d = bwImg(c).*c;
    %fun = @(A,B) ismember(A, B(find(B)));
        
    % create graph
    G = graph(gr, 'upper');
    % perform breadth first search with random start node to get one
    % farthest end
    V1 = bfsearch(G, 1);
    % perform breadth first search with farthest node as start node to get
    % the other farthest end node
    V2 = bfsearch(G, V1(size(V1,1)));
    % create lower sparse matrix for ssp
    grs = sparse(gr');
    % perform shortest path alg. with both farthest end nodes as start and
    % target node in order to get the longest path in DNA fragment. This
    % now is the DNA spine with shorter arms removed
    [dist, path, pred] = graphshortestpath( ... 
        grs, ...
        V1(size(V1,1)), ...
        V2(size(V2,1)), ...
        'Method','Acyclic', ... % we have to ensure that it is acyclic!
        'directed',false);
    singlePath = path;
end


function [newBeginning, newEnd] = elongateDnaBackbone(path, PixelIdxList, bwThick)
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
        endPixel1 = a(path(1));
        endPixel2 = a(path(size(path,2))); 
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
        nextNeighbour1 = mask_end1+endPixel1 == ones(size(mask_end1))*a(path(2));
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
                nextNeighbourIdx = nextNeighbourIdx+nextNeighbourDifference;
                if ( (nextNeighbourIdx < 1) || (nextNeighbourIdx > imgSize(1)*imgSize(2)) )
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
            ones(size(mask_end2))*a(path(size(path,2)-1));
        nextNeighbour2 = fliplr(nextNeighbour2);
        nextNeighbourDifference = mask_end2(find(nextNeighbour2));
        if( nextNeighbourDifference ~= 0)
            nextNeighbourIdx = endPixel2+nextNeighbourDifference;
            while (bwThick(nextNeighbourIdx))
                newEnd = [nextNeighbourIdx; newEnd];
                nextNeighbourIdx = nextNeighbourIdx+nextNeighbourDifference;
                if ( (nextNeighbourIdx < 1) || ...
                        (nextNeighbourIdx > imgSize(1)*imgSize(2)) )
                    break;
                end
            end
        end
    end
end