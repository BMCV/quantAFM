
cpuCores = 4;

running = gcp('nocreate');
if running == 0;
    parpool('local', cpuCores);
end
currentImageDir = 'pictures\p_88A\*.tif';
imageFolderObj = dir(currentImageDir);
imageCount = size(dir(currentImageDir),1);
imageList = cell(1,imageCount);



for index = 1:imageCount
    [image,colorMap] = imread(imageFolderObj(index).name);
    imageList{index} = Image(image, colorMap);
    % Use some kind of threshhold algorithm
    imageList{index}.bwImage = im2bw(imageList{index}.image,0.6);
    imageList{index}.connected =  bwconncomp(imageList{index}.image);
    imageList{index}.region = regionprops(imageList{index}.connected); 
end
