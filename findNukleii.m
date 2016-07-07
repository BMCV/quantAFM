function [centers, radii] = findNukleii(img_bin, img)

warning off
[centers,radii,metric] = imfindcircles(img,[4 6],'ObjectPolarity','bright','Sensitivity',0.96,'EdgeThreshold',0.3);
warning on
if numel(centers)~=0
    centers1 = round(centers);

    % img_bin(diag(img_bin(centers(:,1),centers(:,2)))==0)=0;
    for i=1:size(centers,1)
        if img_bin(centers1(i,2),centers1(i,1)) == 0
            centers(i,:) = 0;
            radii(i) = 0;
        end
    end
    centers(centers(:,1) == 0,:) = [];
    radii(radii==0) =[];
end
end