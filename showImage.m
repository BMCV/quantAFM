function done = showImage(imgObj, imgName, showBB, showThin, purgeInvalid )
 
    fig = figure; set(fig, 'Visible', 'on');
    hold on;
    

    fused = imfuse(imgObj.rawImage, imgObj.bwImgThickDna);
    if (showThin)
        fused = imfuse(fused, imgObj.bwImgThinnedRemoved);
    end
        
    imshow(fused);
    viscircles(imgObj.centers, imgObj.radii, ...
        'LineWidth', 0.5, 'EnhanceVisibility', false);
    rad = ones(length(imgObj.region),1);
    viscircles(imgObj.region, rad,...
        'LineWidth', 0.5, 'EnhanceVisibility', false);

    % for purged files, the purge object is transposed, so take another
    % dimension.
    if (purgeInvalid)
        text(imgObj.region(imgObj.purged,1), imgObj.region(imgObj.purged,2), ...
         num2str ((1: size(imgObj.purged,2))'), 'FontSize', 6);
    else
        text(imgObj.region(:,1), imgObj.region(:,2), ...
         num2str ((1: size(imgObj.region,1))'), 'FontSize', 6);
    end
    
    if (showBB)
        for dnaIndex=1:size(imgObj.dnaList,2)
            curr = imgObj.dnaList{dnaIndex};
            if curr.isValid == 1
                thisBB = curr.bBox.BoundingBox;
                rectangle('Position', thisBB,...
                    'EdgeColor','yellow', 'LineWidth', 1)
            end
        end
    end
    hold off;
    imgName2 = [imgName '_detectedObjects.png'];
    saveas(gcf,imgName2);
    %k = waitforbuttonpress();
    close;
end