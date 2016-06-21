function [angle1,angle2] = measure_angle(dna)
%% 
% dna: dna_bound, dna.attachedNukleo, dna.attachedNukleo.rad, center
% imsize: size of the original image. used to convert linear indices
angle1 = 0;
angle2 = 0;


if numel(dna.attachedNukleo) > 1
    % more than one center!
    return
elseif numel(dna.attachedNukleo) == 0
        return
end

% coordinates of DNA
[x,y] = ind2sub(dna.sizeImg,dna.connectedThinned);
coords = [x,y];


% Distance of each DNA pixel to center of nukleus
D = pdist2(coords,dna.attachedNukleo{1}.localCenter);
D1 = D;

%% case 1: angle between intersecting pixel and centerpoint

    coords(D1<dna.attachedNukleo{1}.rad,:) = [];
    D1(D1<dna.attachedNukleo{1}.rad) = [];
    if size(coords,1) >1
        [~,idx_min] = sort(D1);
        % interesecting points
        p = coords(idx_min(1:2),:);

        angle1 = atan2d(abs(det([p(1,:)-dna.attachedNukleo{1}.localCenter;p(2,:)-dna.attachedNukleo{1}.localCenter])),...
            dot(p(1,:)-dna.attachedNukleo{1}.localCenter,p(2,:)-dna.attachedNukleo{1}.localCenter));
    end

%% case 2: angle between lines fitted on arms
    x1 = x;
    y1 = y;
    % Delete everything inside of nukleus and outside of "outside circle"
    x1(D<dna.attachedNukleo{1}.rad|D>dna.attachedNukleo{1}.rad+5) = [];
    y1(D<dna.attachedNukleo{1}.rad|D>dna.attachedNukleo{1}.rad+5) = [];

    % Create image with two (or unfortunately more) DNA strings to fit line to
    arms = false(dna.sizeImg);
    arms(sub2ind(dna.sizeImg, x1, y1)) = true;

    % workaround more than one arm
    arms_ = bwconncomp(arms);

    if arms_.NumObjects > 2
%         numPixels = cellfun(@numel,arms_.PixelIdxList);
%         [~,idx] = sort(numPixels,'descend');
%         arms( arms_.PixelIdxList{idx(3:end)}) = 0;
%         arms_ = bwconncomp(arms);
    end
    warning off
    if numel(arms_.PixelIdxList) ==2 && size(arms_.PixelIdxList{1},1)>1 && size(arms_.PixelIdxList{2},1)>1
        % get coordinates of arms
        [x1,y1] = ind2sub(dna.sizeImg,arms_.PixelIdxList{1});
        [x2,y2] = ind2sub(dna.sizeImg,arms_.PixelIdxList{2});
        p1 = polyfit(x1,y1,1);
        p2 = polyfit(x2,y2,1);
        if p1(1) == p2(1)
            angle2 = 180;
        else
            p = [x1(1),polyval(p1,x1(1));x2(end), polyval(p2,x2(end))];
            x_intersect = fzero(@(x) polyval(p1-p2,x),3);
            y_intersect = polyval(p1,x_intersect);
            centers = [x_intersect, y_intersect];
            angle2 = atan2d(abs(det([p(1,:)-centers;p(2,:)-centers])),dot(p(1,:)-centers,p(2,:)-centers));            
%             close all
%             imshow(arms);
%             hold on;
%             scatter(p(:,2),p(:,1), 'o','red')
%             plot(polyval(p1,x1),x1);
%             plot(polyval(p2,x2),x2);
%             scatter(centers(2),centers(1),'x','cyan')
%             scatter(dna.attachedNukleo{1}.localCenter(2),dna.attachedNukleo{1}.localCenter(1),'o','blue')
%              disp('---')
%              disp(angle1)
%              disp(angle2)
%             pause
        end
    end
    warning on

end

