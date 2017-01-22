function done = showImage(imgObj, imgName)
 
    %     imshowpair(imageList{index}.rawImage, imageList{index}.bwImgThickDna);
        %imshow(label2rgb(labelmatrix(imageList{index}.connectedThinnedDna));
    fig = figure; set(fig, 'Visible', 'on');
    hold on;
    imshow(imfuse(imgObj.rawImage, imgObj.bwImgThickDna));
    viscircles(imgObj.centers, imgObj.radii, ...
        'LineWidth', 0.5, 'EnhanceVisibility', false);
    rad = ones(length(imgObj.region),1);
    viscircles(imgObj.region, rad,...
        'LineWidth', 0.5, 'EnhanceVisibility', false);
    %     text(imgObj.region(:,1), imgObj.region(:,2),...
    %         [repmat( ' ', size(imgObj.region,1)),num2str((1:1:size(imgObj.region,1))')],'Color', 'b') ;
    % text(imgObj.region(:,1), imgObj.region(:,2),num2str([1:size(imgObj.region,1)]'));
    %text(imgObj.region(:,1), imgObj.region(:,2), num2str ([1: size(imgObj.region,1)]'));

     text(imgObj.region(:,1), imgObj.region(:,2), ...
         num2str ([1: size(imgObj.region,1)]'), 'FontSize', 6);
     %text(imgObj.centers(:,1), imgObj.centers(:,2), num2str(imgObj.attachedDNA));
    for dnaIndex=1:size(imgObj.dnaList,2)
        curr = imgObj.dnaList{dnaIndex};
        if curr.isValid == 1
            thisBB = curr.bBox.BoundingBox;
            rectangle('Position', thisBB,...
                'EdgeColor','yellow', 'LineWidth', 1)
        end
    end
     hold off;
     %imgName = ['../dokumentation/Ergebnisse/alleBilder_auswertungen/' imgName '_detectedObjects.ps'];
      imgName1 = [imgName '_detectedObjects.ps'];
     %print(imgName, '-deps', '-tiff', '-loose', 'driver', '-dpsc');
     if( strcmp(getenv('OS'),'Windows_NT'))
        print(imgName1, '-dwinc', '-r300');
     else
        print(imgName1, '-dpsc', '-r300');
     end
     imgName2 = [imgName '_detectedObjects.png'];
     saveas(gcf,imgName2);
     close;
end