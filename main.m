
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
threshAlgo = 'otsu';
threshAlgo1 = 'intermeans';



parfor index = 1:imageCount
    [image,colorMap] = imread(imageFolderObj(index).name);
%     imageList{index}.metaImage = imfinfo(imageFolderObj(index).name);
    imageList{index} = Image(image, colorMap);
    % Use some kind of threshhold algorithm

    imageList{index}.medImage = medfilt2(imageList{index}.image,[3 3]);
    
    background = imopen(imageList{index}.medImage, strel('disk',15));
    imageList{index}.medImage = imageList{index}.medImage - background;
    thresh = threshold(threshAlgo, imageList{index}.medImage);
    imageList{index}.bwImage = im2bw(imageList{index}.medImage,thresh);
    
    

    fftImage = fftshift(fft2(imageList{index}.medImage));
    imageList{index}.fftImage = mat2gray(log(abs(fftImage)));
    t = threshold(threshAlgo,imageList{index}.fftImage);
    imageList{index}.fftbwImage = im2bw(imageList{index}.fftImage,t);
    [imageList{index}.cArray.center, imageList{index}.cArray.rad] = imfindcircles(imageList{index}.fftImage,500);
    
    
    imageList{index}.connected =  bwconncomp(imageList{index}.bwImage);
    imageList{index}.region = regionprops(imageList{index}.connected);
    
    imageList{index}.bwImageremoved = bwareafilt(imageList{index}.bwImage, [300,900]);
    imageList{index}.bwImageremoved = imageList{index}.bwImage - imageList{index}.bwImageremoved;
    imageList{index}.complementImage = imcomplement(imageList{index}.bwImageremoved);
    imageList{index}.cleanImage = uint8(imageList{index}.bwImageremoved) * uint8(floor(thresh *255*0.8)) + imageList{index}.medImage .* uint8(imageList{index}.complementImage);
    t = threshold(threshAlgo , imageList{index}.cleanImage);
    imageList{index}.bwImage1 = im2bw(imageList{index}.cleanImage,t);
    imageList{index}.bwImageremoved = bwareafilt(imageList{index}.bwImage, [300,900]);
    
end
