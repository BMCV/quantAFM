function [gr, singlePath ]= makeGraph(pxlIdxList, bwImg, bwThick)
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
    % x: don't care, 
    % idx: current pixel index
    % 1: check in bwImgThinnedDna image whether this pixel is 1 or 0
    mask = [1, imgSize(1)-1, imgSize(1), imgSize(1)+1 ];
    % iterate over each entry in pixelIdxList
    for i = 1:size(a,1)
        % apply mask and get indices where mask is true
        res = bwImg(a(i)+mask).*(mask+a(i));
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
    mask2 = [-imgSize(1)-1, -imgSize(1), -imgSize(1)+1, -1, ...
              1, imgSize(1)-1, imgSize(1), imgSize(1)+1 ];
    newBeginning = [];
    newEnd = [];
    if ~(size(a,1) < 2)
        endPixel1 = a(path(1));
        endPixel2 = a(path(size(path,2))); 
        % elongate path with missing pixels at one end
        mask_end1 = mask2;
        mask_end1(endPixel1+mask_end1 < 1 | mask_end1+endPixel1 > imgSize(1)*imgSize(2) ) = [];
        nextNeighbour1 = mask_end1+endPixel1 == ones(size(mask_end1))*a(path(2));
        nextNeighbour1 = fliplr(nextNeighbour1);
        nextNeighbourDifference = mask_end1(find(nextNeighbour1));
        nextNeighbourIdx = endPixel1+nextNeighbourDifference;
    %    if ( ~(nextNeighbourIdx < 1) && ~(nextNeighbourIdx > imgSize(1)*imgSize(2)) )
            while (bwThick(nextNeighbourIdx))
                newBeginning = [nextNeighbourIdx; newBeginning];
                nextNeighbourIdx = nextNeighbourIdx+nextNeighbourDifference;
                if ( (nextNeighbourIdx < 1) | (nextNeighbourIdx > imgSize(1)*imgSize(2)) )
                    break;
                end
            end
    %    end
        mask_end2 = mask2;
        mask_end2(endPixel2+mask_end2 < 1 | mask_end2+endPixel2 > imgSize(1)*imgSize(2) ) = [];
    %    nextNeighbour2 = bwImg(endPixel2+mask_end2);
        nextNeighbour2 = mask_end2+endPixel2 == ones(size(mask_end2))*a(path(size(path,2)-1));
        nextNeighbour2 = fliplr(nextNeighbour2);
        nextNeighbourDifference = mask_end2(find(nextNeighbour2));
        nextNeighbourIdx = endPixel2+nextNeighbourDifference;
        while (bwThick(nextNeighbourIdx))
            newEnd = [nextNeighbourIdx; newEnd];
            nextNeighbourIdx = nextNeighbourIdx+nextNeighbourDifference;
            if ( (nextNeighbourIdx < 1) | (nextNeighbourIdx > imgSize(1)*imgSize(2)) )
                break;
            end
        end
    end
    singlePath =  [ newBeginning; a(path); newEnd];