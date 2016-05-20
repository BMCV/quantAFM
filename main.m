
if( strcmp(getenv('OS'),'Windows_NT'))
    
    addpath(genpath('..\pictures'));
    currentImageDir = '..\pictures\p_Wildtyp\*.tif';
    
else
    addpath(genpath('../pictures'));
    currentImageDir = '../pictures/_Wildtyp/*.tif';
end
        
    
cpuCores = 4;

running = gcp('nocreate');
if running == 0;
    parpool('local', cpuCores);
end
imageFolderObj = dir(currentImageDir);
imageCount = size(dir(currentImageDir),1);
imageList = cell(1,imageCount);



parfor index = 1:imageCount
    [image,colorMap] = imread(imageFolderObj(index).name);
%     imageList{index}.metaImage = imfinfo(imageFolderObj(index).name);
    imageList{index} = Image(image, colorMap);
    % Use some kind of threshhold algorithm
    imageList{index}.medImage = medfilt2(imageList{index}.image,[3 3]);
    threshold = graythresh(imageList{index}.medImage);
    imageList{index}.bwImage = im2bw(imageList{index}.medImage,threshold);
    
    

    fftImage = fftshift(fft2(imageList{index}.image));
    imageList{index}.fftImage = mat2gray(log(abs(fftImage)));
    t = graythresh(imageList{index}.fftImage);
    imageList{index}.fftbwImage = im2bw(imageList{index}.fftImage,t);
    [imageList{index}.cArray.center, imageList{index}.cArray.rad] = imfindcircles(imageList{index}.fftImage,500);
    
    
    imageList{index}.connected =  bwconncomp(imageList{index}.bwImage);
    imageList{index}.region = regionprops(imageList{index}.connected);
    
    imageList{index}.bwImageremoved = bwareafilt(imageList{index}.bwImage, [300,900]);
    imageList{index}.bwImageremoved = imageList{index}.bwImage - imageList{index}.bwImageremoved;
    imageList{index}.complementImage = imcomplement(imageList{index}.bwImageremoved);
    imageList{index}.image = uint8(imageList{index}.bwImageremoved) * uint8(floor(threshold *255*0.8)) + imageList{index}.image .* uint8(imageList{index}.complementImage);
    t = graythresh(imageList{index}.image);
    imageList{index}.bwImage = im2bw(imageList{index}.image,t);
    imageList{index}.bwImageremoved = bwareafilt(imageList{index}.bwImage, [300,900]);
    
end
