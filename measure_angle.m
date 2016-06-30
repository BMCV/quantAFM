function [angle1,angle2] = measure_angle(dna)
%% 
% dna: dna_bound, dna.attachedNukleo, dna.attachedNukleo.rad, center

angle1 = 0;
angle2 = 0;


if numel(dna.attachedNukleo) ~= 1
    % more than one center or no center!
    return
end

% cut out nucleo
dnaImgThinned = dna.bwImageThinned;
mask = bsxfun(@plus, ((1:dna.sizeImg(2)) - dna.attachedNukleo{1}.localCenter(1)).^2, ((transpose(1:dna.sizeImg(1)) - dna.attachedNukleo{1}.localCenter(2)).^2)) < dna.attachedNukleo{1}.rad^2;
dnaImgThinned(mask) = 0;

% Get the two branches
arms = bwconncomp(dnaImgThinned);

% If there are more or less than 2 arms, return
if  numel(arms.PixelIdxList) ~= 2
    return
end
    
%% case 1: angle between intersecting pixel and centerpoint
% get intersecting points
[x1,y1] = ind2sub(dna.sizeImg,arms.PixelIdxList{1});
[x2,y2] = ind2sub(dna.sizeImg,arms.PixelIdxList{2});
[~, idx1] = pdist2([y1,x1], dna.attachedNukleo{1}.localCenter,'euclidean', 'Smallest',1);
[~, idx2] = pdist2([y2,x2], dna.attachedNukleo{1}.localCenter,'euclidean','Smallest',1);
p = [y1(idx1),x1(idx1); y2(idx2), x2(idx2)];
% get lines throug points & nucleo center
a = dna.attachedNukleo{1}.localCenter-p(1,:);
b = dna.attachedNukleo{1}.localCenter- p(2,:);
% get angle
angle1 = acosd(a*b'/(norm(a)*norm(b)));
% show angle
% showAngle


%% case 2: angle between lines fitted on arms
% Delete everything outside of "outside circle"
mask = bsxfun(@plus, ((1:dna.sizeImg(2)) - dna.attachedNukleo{1}.localCenter(1)).^2, ((transpose(1:dna.sizeImg(1)) - dna.attachedNukleo{1}.localCenter(2)).^2)) < (dna.attachedNukleo{1}.rad+5)^2;
dnaImgThinned = dnaImgThinned.*mask;

% Get the two branch stumps
arms = bwconncomp(dnaImgThinned);


if  size(arms.PixelIdxList{1},1)>1 && size(arms.PixelIdxList{2},1)>1
    % get coordinates of arms
    [x1,y1] = ind2sub(dna.sizeImg,arms.PixelIdxList{1});
    [x2,y2] = ind2sub(dna.sizeImg,arms.PixelIdxList{2});
    % Points furthest away from nucleo center
    [~, idx1] = pdist2([y1,x1], dna.attachedNukleo{1}.localCenter,'euclidean', 'Largest',1);
    [~, idx2] = pdist2([y2,x2], dna.attachedNukleo{1}.localCenter,'euclidean','Largest',1);

    % fit lines to branch stumps
    p1 = polyfit(x1,y1,1);
    p2 = polyfit(x2,y2,1);
    if p1(1) == p2(1)
        angle2 = 180;
    else
        p = [polyval(p1,x1(idx1)),x1(idx1); polyval(p2,x2(idx2)),x2(idx2)];
        % get intersection of the two lines
        x_intersect = fzero(@(x) polyval(p1-p2,x),3);
        y_intersect = polyval(p1,x_intersect);
        centers = [y_intersect, x_intersect];
        a = centers-p(1,:);
        b = centers- p(2,:);
        % get angle
        angle2 = acosd(a*b'/(norm(a)*norm(b)));
        % show Angle
        % showAngle2
    end
end


end


