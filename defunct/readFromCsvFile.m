function done = readFromCsvFile(path)
    t146 = xlsread(path);
    len = size(t146,1);
    centers = zeros(len,2);
    centers(:,2) = 1280-t146(:,5);
    centers(:,1) = t146(:,3);
    radii = ones(1,len)*3;
    
    img = imfuse(imageList{23}.rawImage, imageList{23}.bwImgThinnedDna);
    imshow(img);
    hold on;
    viscircles(centers, radii, ...
        'LineWidth', 0.5, 'EnhanceVisibility', false);
    
    text(centers(:,1), centers(:,2), ...
        num2str ([1: len]'), 'FontSize', 6);
    done = 'done'
end