function done = showImage(imgObj, imgName)
 
    %     imshowpair(imageList{index}.rawImage, imageList{index}.bwImgThickDna);
        %imshow(label2rgb(labelmatrix(imageList{index}.connectedThinnedDna));
    fig = figure; set(fig, 'Visible', 'off');
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
     hold off;
     %imgName = ['../dokumentation/Ergebnisse/alleBilder_auswertungen/' imgName '_detectedObjects.ps'];
     imgName = [imgName '_detectedObjects.ps'];
     %print(imgName, '-deps', '-tiff', '-loose', 'driver', '-dpsc');
     if( strcmp(getenv('OS'),'Windows_NT'))
        print(imgName, '-dwinc', '-r300');
     else
        print(imgName, '-dpsc', '-r300');
     end
     close;
end