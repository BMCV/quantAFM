function done = correlateResults(idx, a, b, imgObj)
    
    path1 = ['../dokumentation/Ergebnisse/fastThinning_ChristianLength/denoised_' idx '_fast_ChrLen.csv'];
    path2 = ['../handauswertungen/' idx '.xlsx'];
    target = ['correlation_' idx '.csv'];
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
    [M,I] = min(dists);
    numbers = 1:size(I,2);
    correlation = [I;numbers];
    csvwrite(target, correlation);
    
    showResults(t2, a,b, len, imgObj);
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