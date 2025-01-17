function angle = measure_angle(dna, imsize, case)
%% 
% dna: dna_bound, dna.attachedNukleo, dna.attachedNukleo.rad, center
% imsize: size of the original image. used to convert linear indices

if numel(dna.attachedNukleo) > 2
    % more than one center!
elseif numel(dna.attachedNukleo) == 0
        angle = 0;
        return
end

% coordinates of DNA
[x,y] = ind2sub(imsize,dna.connectedThinned);
coords = [x,y];

% Distance of each DNA pixel to center of nukleus
D = pdist2(coords,dna.attachedNukleo.center);


%% case 1: angle between intersecting pixel and centerpoint
if angle_case == 0
    coords(D<dna.attachedNukleo.rad,:) = [];
    D(D<dna.attachedNukleo.rad) = [];
    [~,idx_min] = sort(D);
    % interesecting points
    p = coords(idx_min(1:2),:);

    angle = atan2d(abs(det([p(1,:)-dna.attachedNukleo.center;p(2,:)-dna.attachedNukleo.center])),dot(p(1,:)-dna.attachedNukleo.center,p(2,:)-dna.attachedNukleo.center));
end

%% case 2: angle between lines fitted on arms
if angle_case == 1
    % Delete everything inside of nukleus and outside of "outside circle"
    x(D<dna.attachedNukleo.rad|D>dna.attachedNukleo.rad+5) = [];
    y(D<dna.attachedNukleo.rad|D>dna.attachedNukleo.rad+5) = [];

    % Create image with two (or unfortunately more) DNA strings to fit line to
    arms = false(imsize);
    arms(sub2ind(imsize, y, x)) = true;

    % workaround more than one arm
    arms_ = bwconncomp(arms);

    if arms_.NumObjects > 2
        numPixels = cellfun(@numel,arms_.PixelIdxList);
        [~,idx] = sort(numPixels,'descend');
        arms( arms_.PixelIdxList{idx(3:end)}) = 0;
        arms_ = bwconncomp(arms);
    end

    % get coordinates of arms
    [x1,y1] = ind2sub(size(img),arms_.PixelIdxList{1});
    [x2,y2] = ind2sub(size(img),arms_.PixelIdxList{2});
    p1 = polyfit(x1,y1,1);
    p2 = polyfit(x2,y2,1);
    p = [polyval(p1,x1(1)), x1(1); polyval(p2,x2(2)),x2(2)];
    x_intersect = fzero(@(x) polyval(p1-p2,x),3);
    y_intersect = polyval(p1,x_intersect);
    centers = [x_intersect, y_intersect];
    angle = atan2d(abs(det([p(1,:)-centers;p(2,:)-centers])),dot(p(1,:)-centers,p(2,:)-centers));

end


