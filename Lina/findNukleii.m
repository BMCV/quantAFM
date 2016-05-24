img = imread('../pictures/p_Wildtyp/150.tif');
img = medianfilter(img);
[centers,radii,metric] = imfindcircles(img,[3 6],'ObjectPolarity','bright','Sensitivity',0.98,'EdgeThreshold',0.3);
figure; imshow(img); viscircles(centers, radii, 'EdgeColor','b');

thresh_low = graythresh(img);
[img_thin,img_bin] = thresh_thin(img,thresh_low);


%[centers_b,radii_b,metric_b] = imfindcircles(img_bin,[3 7],'ObjectPolarity','bright','Sensitivity',0.9,'EdgeThreshold',0.5);
figure; imshow(img_bin); viscircles(centers,radii,'EdgeColor','b');

centers1 = round(centers);

for i=1:size(centers,1)
    if img_bin(centers1(i,2),centers1(i,1)) == 0
        centers(i,:) = 0;
        radii(i) = 0;
    end
end
centers(centers(:,1) == 0,:) = [];
radii(radii==0) =[];

figure; imshow(img_bin); viscircles(centers,radii,'EdgeColor','b');
figure; imshow(img_thin); viscircles(centers,radii,'EdgeColor','b');