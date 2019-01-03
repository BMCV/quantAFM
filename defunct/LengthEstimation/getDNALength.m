function [ dnaObj ] = getDNALength( dnaObj, dnaHasNucleos )
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
% EDIT 26-01-17: Further distinction between errors:
% -1: PixelList was empty
% -2: Constraint not fulfilled (cross-section, or loops)
% -3: Total length smaller/larger than unbound length consraints
% -4: Sum of both arms is smaller/larger than bound length constraints

% Around the picture a quad of 0s is inserted so we dont have to check if
% a point is at the edge of the picture
global MINLENGTH_FREE MAXLENGTH_FREE MINLENGTH_BOUND MAXLENGTH_BOUND RECOVERBACKBONE;
picture = logical(padarray(dnaObj.bwImageThinned, [1 1]));
[height, width] = size(picture);
cc = bwconncomp(picture);
fragmentLen = {0};
% dnaObj.isValid = 1;

% if the image is empty return
if(isempty(cc.PixelIdxList))
    fragmentLen{1} = 0;
    dnaObj.isValid = -1;
    dnaObj.length= fragmentLen;
    return;
end;

% get all pixels that are 1
pixels = cc.PixelIdxList{1};
% create the adjacencyMatrix
cases = createAdjMat(picture, pixels);
adjMat = logical(cases);

% initialize cell struct for fragment length(s)


% remove all "L points". That means points that are connected to exactly to
% two points that are connected to each other thereselves. E.g.
% 0 0 0
% 0 1 1
% 0 1 0
% % Therefore we look at a 2x2 clipping and check if the pixel has only 2
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
    end
else
    dnaObj.isValid = -2;
end


% reverse padding
picture = picture(2:height - 1, 2:width - 1);
% picture is 'logical', might be problematic in any other function?
dnaObj.bwImageThinnedRemoved = picture;
% also recalculate the new PixelIdxList
cc = bwconncomp(picture);
dnaObj.connectedThinnedRemoved = cc.PixelIdxList{1};

% elongation process. Using elongateDnaBackbone function from other length
% determination function
if (RECOVERBACKBONE && dnaObj.isValid)
    
%     imshow(picture)
%     dnaObj.isValid
%     eic = bwmorph(picture, 'endpoints');
%     endpoints = find(eic == 1);
%     firstFirst = endpoints(1);
%     firstFirstIndex = find(pixels == endpoints(1));
%     firstSecond = pixels(find(adjMat(firstFirstIndex, :) == 1));
%     firstLast = endpoints(2);
%     firstLastIndex = find(pixels == endpoints(2));
%     secondLast = pixels(find(adjMat(firstLastIndex, :) == 1));
%     figure()
%     picture([firstFirst, firstSecond, secondLast, firstLast]) = 0;
%     imshow(imfuse(picture, eic))

    [newBeginning, newEnd] = elongateDnaBackbone(dnaObj.connectedThinnedRemoved, dnaObj.bwImage);
%     prev = dnaObj.bwImageThinnedRemoved
    picture(newBeginning) = 1;
    picture(newEnd) = 1;
%     imshow(imfuse(prev, imfuse(dnaObj.bwImage, dnaObj.bwImageThinnedRemoved)));
end

% now with elongation
dnaObj.bwImageThinnedRemoved = picture;
% recalculate again for elongated list
cc = bwconncomp(picture);
dnaObj.connectedThinnedRemoved = cc.PixelIdxList{1};

% no or multiple nucleosome(s), so calculate entire fragment's length
fragmentLen{1} = calcKulpaLength(dnaObj.connectedThinnedRemoved, dnaObj.bwImageThinnedRemoved);
%  fragmentLen{1}
if fragmentLen{1}> MAXLENGTH_FREE || fragmentLen{1} < MINLENGTH_FREE
    dnaObj.isValid = -3;
end
% calculate the fragment lengths, if there is exactly 1 nucleo bound to it
if(dnaHasNucleos && numel(dnaObj.attachedNukleo) >= 1) % dna has 1 nucleosome, so calc length for both arms
    arms = getArmsNucleoIntersection(dnaObj);
    for armIdx = 1:arms.NumObjects % for each arm, ...
        currArm = arms.PixelIdxList{armIdx}; % ... get PixelIdxList ...
        fragmentLen{armIdx+1} = calcKulpaLength(currArm, dnaObj.bwImageThinnedRemoved); % ... and calc its length
    end
    
    % EDIT 26-01-17: There was an error here. Only take the sum of the
    % partial arms, not including the total length as well.
    fragmentLen_ = sum([fragmentLen{2:end}]);

    if (fragmentLen_> MAXLENGTH_BOUND || fragmentLen_ < MINLENGTH_BOUND)
        dnaObj.isValid = -4;
    end
end

% hack because of order of function calls
% ADDED 09-01-17: should we even reset the length for invalid ones?
% if dnaObj.isValid == 0
%     fragmentLen{1} = 0;
% end

% set the struct
dnaObj.length = fragmentLen;
% if (dnaObj.isValid == 1)
%     fragmentLen
% end


end
% function end


% could be implemented as a separate function since it is called in
% determineLengthDNA2 as well.
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