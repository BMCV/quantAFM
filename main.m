
if( strcmp(getenv('OS'),'Windows_NT'))
    addpath(genpath('..\pictures'));    
    currentImageDir = '..\pictures\p_Wildtyp\*.tif';
    
else
    addpath(genpath('../pictures'));
    currentImageDir = '../pictures/p_Wildtyp/*.tif';
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
threshAlgo1 = 'moments';



for index = 1:imageCount
%     if or( strcmp(imageFolderObj(index).name , '.'), strcmp(imageFolderObj(index).name, '..') )
%         continue
%     end
%     [image,colorMap] = imread(strcat(currentImageDir, imageFolderObj(index).name));
%     imageList{index}.metaImage = imfinfo(imageFolderObj(index).name);
    [image,colorMap] = imread(imageFolderObj(index).name);
    imageList{index} = Image(image, colorMap);
    % Use some kind of threshhold algorithm

    imageList{index}.medImage = medfilt2(imageList{index}.image,[3 3]);
    
    imageList{index}.background = imopen(imageList{index}.medImage, strel('disk',15));
    imageList{index}.medImage = imageList{index}.medImage - imageList{index}.background;
    
    imageList{index}.medImage = lowPassFilter(imageList{index}.medImage);
    
    thresh = threshold(threshAlgo, imageList{index}.medImage);
    imageList{index}.bwImage = im2bw(imageList{index}.medImage,thresh);
    
    
    

%     fftImage = fftshift(fft2(imageList{index}.medImage));
%     imageList{index}.fftImage = mat2gray(log(abs(fftImage)));
%     t = threshold(threshAlgo,imageList{index}.fftImage);
%     imageList{index}.fftbwImage = im2bw(imageList{index}.fftImage,t);
%     [imageList{index}.cArray.center, imageList{index}.cArray.rad] = imfindcircles(imageList{index}.fftImage,500);
    
    
    imageList{index}.connected =  bwconncomp(imageList{index}.bwImage);
    imageList{index}.region = regionprops(imageList{index}.connected);
    
    imageList{index}.bwImageremoved = bwareafilt(imageList{index}.bwImage, [150,900]);
    imageList{index}.bwImageremoved = imageList{index}.bwImage - imageList{index}.bwImageremoved;
    imageList{index}.complementImage = imcomplement(imageList{index}.bwImageremoved);
    imageList{index}.cleanImage = uint8(imageList{index}.bwImageremoved) * uint8(floor(thresh *255*0.8)) + uint8(imageList{index}.medImage) .* uint8(imageList{index}.complementImage);
    t = threshold(threshAlgo1 , imageList{index}.cleanImage);
    imageList{index}.bwImage1 = im2bw(imageList{index}.cleanImage,t);
    imageList{index}.bwImageremoved = bwareafilt(imageList{index}.bwImage, [150,900]);
    imageList{index}.thinImage = bwmorph(imageList{index}.bwImageremoved,'thin',Inf);
    [centers, radii] = findNukleii(imageList{index}.bwImageremoved, imageList{index}.medImage);
    imageList{index}.centers = centers;
    imageList{index}.radii = radii;
    
    if( strcmp(getenv('OS'),'Windows_NT'))
        
    imwrite(imageList{index}.medImage , ['..\pictures\medImage\' 'med' imageFolderObj(index).name ]);
    imwrite(imageList{index}.complementImage , ['..\pictures\complementImage\'  'com' imageFolderObj(index).name ]);
    imwrite(imageList{index}.cleanImage , ['..\pictures\cleanImage\' 'clean' imageFolderObj(index).name ]);
    imwrite(imageList{index}.bwImageremoved , ['..\pictures\bwImageremoved\' 'bwremoved' imageFolderObj(index).name ]);
    imwrite(imageList{index}.bwImage1 , ['..\pictures\bwImage1\' 'bw1' imageFolderObj(index).name ]);
    imwrite(imageList{index}.bwImage , ['..\pictures\bwImage\' 'bw' imageFolderObj(index).name ]);
    imwrite(imageList{index}.background , ['..\pictures\background\' 'bckground' imageFolderObj(index).name ]);
%     imwrite(imageList{index}.fftImage , ['..\pictures\fftImage\' 'fftImage' imageFolderObj(index).name ]);
    imwrite(imageList{index}.thinImage , ['..\pictures\thinImage\' 'thin' imageFolderObj(index).name ]);
    
%     cd('..\biomedizinischebildanalyse') ;   
%         addpath(genpath('..\pictures\'));

    
    else
    end
    
    
end
