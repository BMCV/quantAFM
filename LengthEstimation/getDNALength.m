function [ result, picture ] = getDNALength( picture )
% @author: Dennis Aumiller, Philip Hausner
% This function evaluates the length of an DNA object in a picture clipping
% The function needs a picture in which exactly one DNA string is visible
% (and no more than one object)
% The return values are the result (so the length of the DNA) and the
% processed picture. The returned picture is the original picture after
% pruning and only the "main string" of the DNA is visible anymore

% Around the picture a quad of 0s is inserted so we dont have to check if
% a point is at the edge of the picture
picture = logical(padarray(picture, [1 1]));
[height, width] = size(picture);
cc = bwconncomp(picture);
% get all pixels that are 1
pixels = cc.PixelIdxList{1};

% create the adjacencyMatrix
cases = createAdjMat(picture, pixels);
adjMat = logical(cases);

% remove all "L points". That means points that are connected to exactly to
% two points that are connected to each other thereselves. E.g.
% 0 0 0
% 0 1 1
% 0 1 0
% Therefore we look at a 2x2 clipping and check if the pixel has only 2
% connections to other pixels

for i = 1:length(pixels)
    % do not check points that are connected to more than two points or to
    % only 1 (endpoint)
    if (sum(adjMat(i,:)) ~= 2 && sum(adjMat(1,:)) ~= 4)
        continue;
    end
    % take 2x2 clipping out of picture
    % 1 4 7
    % 2 5 8
    % 3 6 9
    % this are the values at the respective clipping points as stated
    % above
    a1 = picture(pixels(i) - height - 1);
    a2 = picture(pixels(i) - height);
    a3 = picture(pixels(i) - height + 1);
    a4 = picture(pixels(i) - 1);
    a5 = picture(pixels(i));
    a6 = picture(pixels(i) + 1);
    a7 = picture(pixels(i) + height - 1);
    a8 = picture(pixels(i) + height);
    a9 = picture(pixels(i) + height + 1);
    % create the 2x2 mask
    fragment = [a1,a4;a2,a5];
    if isequal(fragment, [0,1;1,1])
        % delete L point
        picture(pixels(i)) = 0;
        % indicate in pixel list that there was an L point without losing
        % the information about the position
        pixels(i) = -pixels(i);
        % correct adjMat
        adjMat(i,:) = 0;
        adjMat(:,i) = 0;
        cases(i,:) = 0;
        cases(:,i) = 0;
    end
    
    % as above
    fragment = [a2,a5;a3,a6];
    if isequal(fragment, [1,1;0,1])
        picture(pixels(i)) = 0;
        pixels(i) = -pixels(i);
        adjMat(i,:) = 0;
        adjMat(:,i) = 0;
        cases(i,:) = 0;
        cases(:,i) = 0;
    end
    
    % as above
    fragment = [a4,a7;a5,a8];
    if isequal(fragment, [1,0;1,1])
        picture(pixels(i)) = 0;
        pixels(i) = -pixels(i);
        adjMat(i,:) = 0;
        adjMat(:,i) = 0;
        cases(i,:) = 0;
        cases(:,i) = 0;
    end
    
    % as above
    fragment = [a5,a8;a6,a9];
    if isequal(fragment, [1,1;1,0])
        picture(pixels(i)) = 0;
        pixels(i) = -pixels(i);
        adjMat(i,:) = 0;
        adjMat(:,i) = 0;
        cases(i,:) = 0;
        cases(:,i) = 0;
    end
end
clear a1 a2 a3 a4 a5 a6 a7 a8 a9 fragment

% get branch- and endpoints
sic = bwmorph(picture, 'branchpoints');
eic = bwmorph(picture, 'endpoints');
% in branches we save the indices of the branchpoints
branches = find(sic);

% constraint to differ between objects that have a circle in them and
% objects we actually can use to determine the length
constraint = 2;
for i = 1:size(sum(sic(:)))
    if sum(adjMat(find(pixels == branches(1)), :)) == 4
        constraint = constraint + 1;
    end
end

% only evaluate object that are valid. With the constraint we can differ
% between valid objects and objects that have circles in them
if (sum(eic(:)) - sum(sic(:))) >= constraint
    % now we want to delete branches out of our graph until we only have our
    % longest branch left
    % Therefore we first calculate the distance from every branchpoint to every
    % other point and delete the way to the nearest endpoint (or from that
    % endpoint to the next branchpoint it reaches)
    % iterate while we have branchpoints remaining
    while(~isempty(branches))
        % create a mask in which only the branchpoint we want to look at is
        % marked
        singleBranchPointMask = false(size(picture));
        singleBranchPointMask(branches(1)) = 1;
        % estimate distances to all other points with bwdistgeodesic, we
        % just need the distance to endpoints later though
        distances = bwdistgeodesic(picture, singleBranchPointMask);
        % here every distance to a point that is not part of our object is NaN,
        % so we have to check the indices for the endpoints later on with
        % find(distances > 0)
        distances = distances .* eic;
        % get nearest endpoint
        minDistance = min(distances(find(distances > 0)));
        % nearestEndpoints is necessary if there are two (or more) 
        % endpoints that have the same distance to the branchpoint
        nearestEndpoints = find(distances == minDistance);
        endPointIndex = find(pixels == nearestEndpoints(1));
        % now delete all pixels to the next branchpoint (not necessarily
        % the branchpoint we started to look at!
        % For that it has to be deleted out of the picture and the adjMat
        % has to be altered
        % also the endPoint and the branchpoint we reach get deleted (in
        % the masks eic and sic and in the branches list
        % endPointIndex is always the index of the point that would
        % currently be an endPoint and not of the starting endpoint
        eic(pixels(endPointIndex)) = 0;
        while (sic(pixels(endPointIndex)) == 0) 
            % delete pixel in picture
            picture(pixels(endPointIndex)) = 0;
            % if only one point can be reached from the current point
            if (sum(adjMat(endPointIndex,:)) == 1)
                % get index of current endPoint and save it temporarily
                % since we later have to go to the next point but need to
                % remember which endPoint has to be deleted (which has to
                % be done after going to the next point)
                tempIndex = find(adjMat(endPointIndex,:));
            % more than one point can be reached
            else
                % all indices that can be reached
                indicesTemp = find(adjMat(endPointIndex,:));
                % iterate over all reachable points and take the one that
                % is a branchpoint. If we would take another of the pixels
                % the for loop wouldnt break and we would destroy the
                % connected DNA strand
                for j = 1:length(indicesTemp)
                    % detect if branchpoint
                    if sic(pixels(indicesTemp(j))) == 1
                        % set the index to the index of the branchpoint
                        tempIndex = indicesTemp(j);
                        break;
                    end
                end
                clear indicesTemp
            end
            % alter adjMat so the deleted point cant be reached anymore
            adjMat(endPointIndex,:) = 0;
            adjMat(:,endPointIndex) = 0;
            cases(endPointIndex,:) = 0;
            cases(:,endPointIndex) = 0;
            % negate the value of the endPoint in the pixellist so we
            % indicate the point is deleted without losing information
            pixels(endPointIndex) = -pixels(endPointIndex);
            % set new endpoint to iterate over
            endPointIndex = tempIndex;
        end
        % delete branchpoint out of the branchpointindexlist branches and
        % out of the matrix sic in which the place of the branchpoint is
        % saved
        branches = branches(branches ~= branches(find(branches == pixels(endPointIndex))));
        sic(pixels(endPointIndex)) = 0;
    end
end
clear i j tempIndex branches sic singleBranchPointMask minDistance

% now calculate length via Kulpa estimation
% sum up number of pixels that are connected diagonally
odd = sum(sum(mod(cases,2))) / 2;

% even is the number of pixels that are connected horizontically or
% vertically
even = (sum(sum(adjMat)) / 2) - odd;

% evaluate the length using the Kulpa estimator
result = 0.948 * even + 1.343 * odd;

% reverse padding
picture = picture(2:height - 1, 2:width - 1);

end

