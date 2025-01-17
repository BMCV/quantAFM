function [angle1,angle2] = measure_angle(dna)
%% 
% dna: dna_bound, dna.attachedNukleo, dna.attachedNukleo.rad, center
warning off;
angle1 = 0;
angle2 = 0;

global ANGLERADIUS ANGLETHRESHOLD GROWING;

    % more than one center or no center 
if numel(dna.attachedNukleo) ~= 1    
    return
end
    % no DNA backbone present!
if size(dna.connectedThinnedRemoved,1) < 2
    dna.isValid = 0;
    return
end

% cut out nucleo
dnaImgThinned = dna.bwImageThinnedRemoved;
mask = bsxfun(@plus, ((1:dna.sizeImg(2)) - dna.attachedNukleo{1}.localCenter(1)).^2, ...
             ((transpose(1:dna.sizeImg(1)) - dna.attachedNukleo{1}.localCenter(2)).^2)) ...
             < dna.attachedNukleo{1}.rad^2;
dnaImgThinned(mask) = 0;

% Get the two branches
arms = bwconncomp(dnaImgThinned);

% If there are more or less than 2 arms, return
if  numel(arms.PixelIdxList) == 2

    
%% case 1: angle between intersecting pixel and centerpoint
% get intersecting points
[x1,y1] = ind2sub(dna.sizeImg,arms.PixelIdxList{1});
[x2,y2] = ind2sub(dna.sizeImg,arms.PixelIdxList{2});
[~, idx1] = pdist2([y1,x1], dna.attachedNukleo{1}.localCenter,'euclidean', 'Smallest',1);
[~, idx2] = pdist2([y2,x2], dna.attachedNukleo{1}.localCenter,'euclidean','Smallest',1);
p = [y1(idx1),x1(idx1); y2(idx2), x2(idx2)];
% get lines through points & nucleo center
a = dna.attachedNukleo{1}.localCenter-p(1,:);
b = dna.attachedNukleo{1}.localCenter- p(2,:);
% get angle
angle1 = real(acosd(a*b'/(norm(a)*norm(b))));
% show angle
% showAngle

end

%% case 2: angle between lines fitted on arms
% Delete everything outside of "outside circle"


% need this for "growing"
dnaImgThinned_org = dnaImgThinned;
radius = dna.attachedNukleo{1}.rad+ ANGLERADIUS;
mask = bsxfun(@plus, ((1:dna.sizeImg(2)) - dna.attachedNukleo{1}.localCenter(1)).^2, ...
            ((transpose(1:dna.sizeImg(1)) - dna.attachedNukleo{1}.localCenter(2)).^2)) ...
            <radius^2;
dnaImgThinned = dnaImgThinned.*mask;
arms = bwconncomp(dnaImgThinned);

% When "growing" is set, increase the radius of the circle until we
% intersect the object.
if (GROWING)
    % Start with 1.5x the original radius enlargement.
    i = 1;
    while numel(arms.PixelIdxList) < 2 && i<4
        % We now have two separate masks. The inner one "grows" the
        % nucleosome radius, whereas the outer one determines the elements
        % to include when looking for the elements that are used to
        % determine fitting.
        radius = i*ANGLERADIUS/2 + dna.attachedNukleo{1}.rad;
        mask1 = bsxfun(@plus, ((1:dna.sizeImg(2)) - dna.attachedNukleo{1}.localCenter(1)).^2, ...
                ((transpose(1:dna.sizeImg(1)) - dna.attachedNukleo{1}.localCenter(2)).^2)) ...
                <radius^2;
        mask2 = bsxfun(@plus, ((1:dna.sizeImg(2)) - dna.attachedNukleo{1}.localCenter(1)).^2, ...
                ((transpose(1:dna.sizeImg(1)) - dna.attachedNukleo{1}.localCenter(2)).^2)) ...
                <(radius+ANGLERADIUS)^2;
        dnaImgThinned = dnaImgThinned_org;
        dnaImgThinned(mask1) = 0;
        dnaImgThinned = dnaImgThinned.*mask2;
        arms = bwconncomp(dnaImgThinned);
        i = i+1;
        
%         showAngle3(dna, mask1, mask2, radius);
    end
end

% Return in both cases if we still have a wrong number of elements.
if  numel(arms.PixelIdxList) ~= 2
    return
end

if  size(arms.PixelIdxList{1},1)>1 && size(arms.PixelIdxList{2},1)>1
    % get coordinates of arms
    [x1,y1] = ind2sub(dna.sizeImg,arms.PixelIdxList{1});
    [x2,y2] = ind2sub(dna.sizeImg,arms.PixelIdxList{2});
    % Points furthest away from nucleo center
    [~, idx1] = pdist2([y1,x1], dna.attachedNukleo{1}.localCenter,'euclidean', 'Largest',1);
    [~, idx2] = pdist2([y2,x2], dna.attachedNukleo{1}.localCenter,'euclidean','Largest',1);
    
    % TODO: differentiate between different cases. Fit line to whatever dimension
    % has more variance. Otherwise, the fit will be extremely bad!
    p1 = polyfit(x1,y1,1);
    p2 = polyfit(x2,y2,1);

    
    % now a bit more complicated
    if p1(1) == p2(1)
        angle2 = 180;
    elseif isinf(p1(1))
        p1(1) = sign(p1(1))*1e10;
    elseif isinf(p2(1))
        p2(1) = sign(p2(1))*1e10;
    elseif (isnan(p1(1)) || isnan(p2(1)) )
        angle2 = 0;
    else
        p = [polyval(p1,x1(idx1)),x1(idx1); polyval(p2,x2(idx2)),x2(idx2)];
        % get intersection of the two lines
        x_intersect = fzero(@(x) polyval(p1-p2,x),3);
        y_intersect = polyval(p1,x_intersect);
        centers = [y_intersect, x_intersect];
        a = centers-p(1,:);
        b = centers-p(2,:);
        % get angle
        angle2 = acosd(a*b'/(norm(a)*norm(b)));
        % Under certain circumstances, small angles can be "inverted", i.e.
        % just 180-X the degree we actually want. In that case, we simply
        % invert it again. A rule of thumb is that we can process it in a
        % way were the polyfit is pointing towards a direction that is
        % "somewhat similar" (i.e. dotproduct > 0).
        if angle2 < ANGLETHRESHOLD
            % show Angle
            angle_dot = dot(p1/sum(p1), p2/sum(p2));
            if angle_dot < 0
                angle2 = 180 - angle2;
            end
%         showAngle2(dna, mask, p, centers, p1, p2, x1, y1, x2, y2, flip1, flip2)     
        end
    end
end
warning on;

end


