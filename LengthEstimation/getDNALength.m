function [ dnaObj ] = getDNALength( dnaObj )
% @author: Dennis Aumiller, Philip Hausner
% This function evaluates the length of an DNA object in a picture clipping
% The function needs a picture in which exactly one DNA string is visible
% (and no more than one object)
% UPDATE: NOW THE FUNCTION TAKES A DNA OBJECT WHICH HAS TO CONTAIN THIS
% THINNED DNA IMAGE
% The return value is the object itself which now has his length attribute
% updated and -1 if the length could not be determined. Secondly,
% bwImageThinnedRemoved is set to an image that was pruned and has its
% L-points removed.
% Have fun using this function!

% Around the picture a quad of 0s is inserted so we dont have to check if
% a point is at the edge of the picture
picture = logical(padarray(dnaObj.bwImageThinned, [1 1]));
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
    % do not check points that are connected to more than two points
    if (sum(adjMat(i,:)) ~= 2)
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
    if (a2 + a4 + a6 + a8 == 4 || a1 + a3 + a7 + a9 == 4)
        continue
    end
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

% if two strands overlap the constraint is increased --> wont be processed
for i = 1:sum(sic(:))
    if sum(adjMat(find(pixels == branches(i)), :)) == 4
        constraint = inf;
    end
end

% bool to check if the DNA is processed
correctDNA = 1;

% only evaluate object that are valid. With the constraint we can differ
% between valid objects and objects that have circles in them
if (sum(eic(:)) - sum(sic(:))) == constraint
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
            end
            % alter adjMat so the deleted point cant be reached anymore
            adjMat(endPointIndex,:) = 0;
            adjMat(:,endPointIndex) = 0;
            cases(endPointIndex,:) = 0;
            cases(:,endPointIndex) = 0;
            % negate the value of the endPoint in the pixellist so we
            % indicate the point is deleted without losing information
            pixels(endPointIndex) = -pixels(endPointIndex);
            %pixels(endPointIndex) = 0;
            % set new endpoint to iterate over
            endPointIndex = tempIndex;
        end
        
        % delete branchpoint out of the branchpointindexlist branches and
        % out of the matrix sic in which the place of the branchpoint is
        % saved
        branches = branches(branches ~= branches(find(branches == pixels(endPointIndex))));
        sic(pixels(endPointIndex)) = 0;
        % check if this (not anymore) branchpoint is now an L-point
        % get the index in 2 dimensions for the point
        [p1, p2] = ind2sub(size(picture), pixels(endPointIndex));
            % create 3x3 picture
            if (isLPoint(picture(p1 - 1 : p1 + 1, p2 - 1 : p2 + 1)))
                % if it is delete it
                % alter adjMat so the deleted point cant be reached anymore
                adjMat(endPointIndex,:) = 0;
                adjMat(:,endPointIndex) = 0;
                cases(endPointIndex,:) = 0;
                cases(:,endPointIndex) = 0;
                % negate the value of the endPoint in the pixellist so we
                % indicate the point is deleted without losing information
                picture(pixels(endPointIndex)) = 0;
                %pixels(endPointIndex) = -pixels(endPointIndex);
            end
        clear p1 p2
    end
else
    correctDNA = 0;
end
clear i j tempIndex branches sic singleBranchPointMask minDistance

% calculate length if DNA is a valid string
if correctDNA
    % #lifehack
    % retrace the ends of the thinned image to the original image so no
    % length gets lost
    % get the two remaining endpoints and their indices in the image
    endpointzz = find(eic);
    endpointzzIndex1 = find(pixels == endpointzz(1));
    endpointzzIndex2 = find(pixels == endpointzz(2));
    % get the neighbours of the endpointzz
    connectedthingyIndex1 = find(adjMat(:,endpointzzIndex1) == 1);
    connectedthingyIndex2 = find(adjMat(:,endpointzzIndex2) == 1);
    % fuse the four points in one array. The order is important for the
    % elongateDNAbackbone function
    crippledPixels = [endpointzz(1),pixels(connectedthingyIndex1),pixels(connectedthingyIndex2),endpointzz(2)];
    % elongateDnaBackbone retraces the points (thats the central point of
    % this step)
    [newBeginning, newEnd] = elongateDnaBackbone(crippledPixels', logical(padarray(dnaObj.bwImage, [1 1])));
    
    % for length detection these two variables have to be added to the
    % other odd and even values of the rest of the DNA string (see Kulpa
    % for reference)
    OddlySpecial = 0;
    EvenSpecial = 0;
    
    % newBeginning pixels odd or even?
    for index = 2:length(newBeginning)
        testsum = abs(newBeginning(index) - newBeginning(index-1));
        if testsum == 1 || testsum == height
            EvenSpecial = EvenSpecial + 1;
        else
            OddlySpecial = OddlySpecial + 1;
        end
    end
    % newBeginning to Endpoint1 odd or even?
    if length(newBeginning) >= 1
        testsum = abs(endpointzz(1) - newBeginning(end));
        if testsum == 1 || testsum == height
            EvenSpecial = EvenSpecial + 1;
        else
            OddlySpecial = OddlySpecial + 1;
        end
    end
       
    % do the same for the new Endpoints
    for index = 2:length(newEnd)
        testsum = abs(newEnd(index) - newEnd(index-1));
        if testsum == 1 || testsum == height
            EvenSpecial = EvenSpecial + 1;
        else
            OddlySpecial = OddlySpecial + 1;
        end
    end
    
    if length(newEnd) >= 1
        testsum = abs(endpointzz(2) - newEnd(end));
        if testsum == 1 || testsum == height
            EvenSpecial = EvenSpecial + 1;
        else
            OddlySpecial = OddlySpecial + 1;
        end
    end
    
    % now calculate length via Kulpa estimation
    % sum up number of pixels that are connected diagonally
    odd = sum(sum(mod(cases,2))) / 2 + OddlySpecial;

    % even is the number of pixels that are connected horizontically or
    % vertically
    even = (sum(sum(adjMat)) / 2) - odd + EvenSpecial;

    % evaluate the length using the Kulpa estimator
    dnaObj.length = 0.948 * even + 1.343 * odd;
else
    dnaObj.length = -1;
end

% reverse padding
picture = picture(2:height - 1, 2:width - 1);
% picture is 'logical', might be problematic in any other function?
dnaObj.bwImageThinnedRemoved = picture;

end

