function [centers, radii] = findNukleii(img_bin, img)
global MINRADIUS MAXRADIUS SENSITIVITY EDGETHRESHOLD;

warning off
[centers,radii,metric] = imfindcircles(img,[MINRADIUS MAXRADIUS],...
                            'ObjectPolarity','bright',...
                            'Sensitivity',SENSITIVITY,...
                            'EdgeThreshold',EDGETHRESHOLD);
warning on
if numel(centers)~=0
    centers1 = round(centers);

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