function [dnaObj] = determineDnaLength2(dnaObj, dnaHasNucleos)
    currPxlList = dnaObj.connectedThinned;
    bwImgThin = dnaObj.bwImageThinned;
    bwImgThick = dnaObj.bwImage;
    bwImgThinnedRemoved = zeros(size(bwImgThin));
    fragmentLen = {};
    % if fragment too small or too large don't compute the rest
    if( ~(size(currPxlList,1) < 3) ) 
        % create graph from its PixelIdxList and, from that, get the
        % fragment backbone
        [gr, singlePath] = getDnaBackbone(currPxlList, bwImgThin);
        % both end parts of the fragment were lost during the
        % thinning step, so get them back from the thickDNA fragment
        [newBeginning, newEnd] = ... 
            elongateDnaBackbone(currPxlList(singlePath), bwImgThick);
        % ... and add the respective pixels to the DNA backbone
        singlePath = [newBeginning; currPxlList(singlePath); newEnd];
        % create bw image of elongated DNA backbone
        bwImgThinnedRemoved(singlePath) = 1;
        dnaObj.bwImageThinnedRemoved = bwImgThinnedRemoved;
        dnaObj.connectedThinnedRemoved = singlePath;
        
        %% for each arm or for entire fragment, calculate length with Kulpa Estimator
        if(dnaHasNucleos && numel(dnaObj.attachedNukleo) == 1) % dna has 1 nucleosome, so calc length for both arms
            arms = getArmsNucleoIntersection(dnaObj);
            for armIdx = 1:arms.NumObjects % for each arm, ...
                currArm = arms.PixelIdxList{armIdx}; % ... get PixelIdxList ...
                fragmentLen{armIdx} = calcKulpaLength(currArm, bwImgThinnedRemoved); % ... and calc its length
            end
        else % no or multiple nucleosome(s), so calculate entire fragment's length
            fragmentLen{1} = calcKulpaLength(singlePath, bwImgThinnedRemoved);
        end
        dnaObj.connectedThinnedRemoved = singlePath;
        
        % set invalid flag if DNA object does not fit length criteria
        if ((size(singlePath, 1) < dnaObj.MIN_LENGTH) || ...
            (size(singlePath, 1) > dnaObj.MAX_LENGTH))
            dnaObj.isValid = 0;
        else
            dnaObj.isValid = 1;
        end
    else
        fragmentLen{1} = 0;
        dnaObj.isValid = 0;
    end
    dnaObj.length = fragmentLen;
%    imshow(imfuse(bwImgThick, bwImgThinnedRemoved));
%    imwrite(imfuse(bwImgThick, bwImgThinnedRemoved) , '../pictures/DNA_spine_thickDna_overlay.tif');
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
        gr(i,:) = ismember(a, res(res ~=0 ));
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

function fragmentLen = calcKulpaLength(pxlIdxList, bwImgThinnedRemoved)
    % Length = 0.948*Ne + 1.343*No
    % with:     Ne - number of even pixels
    %           No - number of odd pixels
    % Here, even pixels are those that are, in their 8-neighourhood
    % connected along the even axes, and odd pixels are those connected
    % along the uneven axes.
    %   3   2   1        0   1   0      would have:
    %   4  pxl  0   =>   0   1   0  =>  1 odd pixel
    %   5   6   7        1   0   0      1 even pixel
    [row, col]= ind2sub(size(bwImgThinnedRemoved),pxlIdxList');
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
    fragmentLen = 0.948*numberOfEvenPixels + 1.343 * numberOfOddPixels;
end
