function done = correlateResults(idx, a, b, imgObj)
    
    path1 = ['../dokumentation/Ergebnisse/auswertungen/fastThinning_ChristianLength/denoised_' idx '.tif_fast_Chr.csv'];
    path2 = ['../handauswertungen/' idx '.xlsx'];
    target1 = ['correlation1_' idx '.csv'];
    target2 = ['correlation2_' idx '.csv'];
    len = imgObj.imgSize(1);
    t1 = csvread(path1,1,1);
    t2 = xlsread(path2);

    x1 = t1(:,1);
    y1 = t1(:,2);
    x2 = len-t2(:,a);
    y2 = t2(:,b);
    vec1 = [x1, y1];
    vec2 = [x2, y2];
    dists = pdist2(vec1, vec2);
    dists2 = pdist2(vec2, vec1);
    [M1,I1] = min(dists);
    [M2,I2] = min(dists2);
    numbers1 = 1:size(I1,2);
    numbers2 = 1:size(I2,2);
    correlation1 = [I1;numbers1];
    correlation2 = [I2;numbers2];
    csvwrite(target1, correlation1);
    csvwrite(target2, correlation2);
    
    %showResults(t2, a,b, len, imgObj);
end

function showResults(t2, a,b, imgLen, imgObj)
    len = size(t2,1);
    centers = zeros(len,2);
    centers(:,2) = imgLen-t2(:,a);
    centers(:,1) = t2(:,b);
    radii = ones(1,len)*3;
    
    img = imfuse(imgObj.rawImage, imgObj.bwImgThickDna);
    imshow(img);
    hold on;
    viscircles(centers, radii, ...
        'LineWidth', 0.5, 'EnhanceVisibility', false);
    
    text(centers(:,1), centers(:,2), ...
        num2str ([1: len]'), 'FontSize', 6);
    hold off;
end