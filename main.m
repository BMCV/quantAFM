
cpuCores = 4;

running = gcp('nocreate');
if running == 0;
    parpool('local', cpuCores);
end
currentImageDir = 'pictures\p_Wildtyp\*.tif';
imageFolderObj = dir(currentImageDir);
imageCount = size(dir(currentImageDir),1);
imageList = cell(1,imageCount);



parfor index = 1:imageCount
    [image,colorMap] = imread(imageFolderObj(index).name);
%     imageList{index}.metaImage = imfinfo(imageFolderObj(index).name);
    imageList{index} = Image(image, colorMap);
    % Use some kind of threshhold algorithm
    imageList{index}.medImage = medfilt2(imageList{index}.image,[3 3]);
    threshhold = graythresh(imageList{index}.medImage);
    imageList{index}.bwImage = im2bw(imageList{index}.medImage,threshhold);
    
    

    fftImage = fftshift(fft2(imageList{index}.image));
    imageList{index}.fftImage = mat2gray(log(abs(fftImage)));
    t = graythresh(imageList{index}.fftImage);
    imageList{index}.fftbwImage = im2bw(imageList{index}.fftImage,t);
    [imageList{index}.cArray.center, imageList{index}.cArray.rad] = imfindcircles(imageList{index}.fftImage,500);
    
    
    imageList{index}.connected =  bwconncomp(imageList{index}.bwImage);
    imageList{index}.region = regionprops(imageList{index}.connected);
%     imageList{index}.bwImageremoved = bwareaopen(imageList{index}.bwImage, 300);
    imageList{index}.bwImageremoved = bwareafilt(imageList{index}.bwImage, [300,900]);
    
    
    
end
