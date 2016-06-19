function [gr, singlePath ]= makeGraph(pxlIdxList, bwImg)
    a = pxlIdxList;
    imgSize = size(bwImg);
%    gr = diag(ones(1,size(a,1)));
    gr = zeros(size(a,1),size(a,1));
    mask = [1, imgSize(1)-1, imgSize(1), imgSize(1)+1 ];
    for i = 1:size(a,1)

        res = bwImg(a(i)+mask).*(mask+a(i));
        % we need to exclude border cases!
        % so, we test whether the found indices really occur in a
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
    singlePath = a(path);