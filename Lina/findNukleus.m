img = imread('../../pictures/p_Wildtyp/150.tif');
img = medianfilter(img);
[centers,radii,metric] = imfindcircles(img,[3 7],'ObjectPolarity','bright','Sensitivity',0.98,'EdgeThreshold',0.3);
figure; imshow(img); viscircles(centers, radii, 'EdgeColor','b');