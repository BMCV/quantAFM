function angle = measure_angle(img_thin,centers,radii, angle_case)
%% 
% img_thin: cropped image of one DNA
% centers, radii: Centerpoint and radius of nukleus
% angle_case: 0 for angle between arms on nukleus center, 1 for
% angle between fitted lines on arms
% return angle: caluclated angle

if numel(centers) > 2
    % more than one center!
elseif numel(centers) == 0
        angle = 0;
        return
end
% coordinates of DNA
[y,x] = find(img_thin);
coords = [x,y];
% Distance of each DNA pixel to center of nukleus
D = pdist2(coords,centers);


%% case 1: angle between intersecting pixel and centerpoint
if angle_case == 1
    coords(D<radii,:) = [];
    D(D<radii) = [];
    [~,idx_min] = sort(D);
    % interesecting points
    p = coords(idx_min(1:2),:);

    angle = atan2d(abs(det([p(1,:)-centers;p(2,:)-centers])),dot(p(1,:)-centers,p(2,:)-centers));
end

%% case 2: angle between lines fitted on arms
if angle_case == 2
    % Delete everything inside of nukleus and outside of "outside circle"
    x(D<radii|D>radii+5) = [];
    y(D<radii|D>radii+5) = [];

    % Create image with two (or unfortunately more) DNA strings to fit line to
    arms = false(size(img));
    arms(sub2ind(size(img), y, x)) = true;

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


