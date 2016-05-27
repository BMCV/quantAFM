function [ image ] = lowPassFilter( image )
%FFT Summary of this function goes here
%   Detailed explanation goes here

 	image = fftshift(fft2(image));
    gray_image = mat2gray(log(abs(image)));
    threshold = graythresh(gray_image);
    bwImage = im2bw(gray_image,threshold);
    
    SE = strel('disk', 4, 8);
    bwImage(floor(size(bwImage, 1)/2) - 20 : floor(size(bwImage, 1)/2) + 20, ...
    floor(size(bwImage, 2)/2) - 20 : floor(size(bwImage, 2)/2) + 20) = 1;
    bwImage = imerode(bwImage, SE);
    bwImage = bwareafilt(bwImage, 1, 'largest');
    bwImage = imfill(bwImage, 'holes');
    % delete the pixels in every row that has less than 75 white pixels
%     image(sum(image) < 75, :) = 0;
%    for i = 1:size(bwImage,1)
%        if(sum(bwImage(i,:)) < 75)
%            bwImage(i,:) = 0;
%        end
%    end
    % faster:
    % sum(img, 2) sums rows, returns [1,n] array
    % ~(sum(img,2) >= 75) returns inverted binary mask of where that array
    %   is larger than 75
    % find(a) returns indices of non-zero elements
    bwImage(find(~(sum(bwImage, 2)>=75)), :) = 0;
    
    SE = strel('disk', 80, 8);
    bwImage = imdilate(bwImage, SE); 
    cc = regionprops(bwImage, 'Centroid', 'MajorAxisLength', 'MinorAxisLength');
%    imshow(bwImage);
%    hold on;
%    viscircles(cc.Centroid, rad); 
    xcoord = (1:size(bwImage, 2)) - cc.Centroid(1);
    ycoord = (1:size(bwImage, 1)) - cc.Centroid(2);
    
    [y,x] = ndgrid(ycoord, xcoord);
    myZeros = zeros(size(bwImage));
    
    if (size(bwImage,1) == size(bwImage,2))
        rad = (0.55 * cc.MajorAxisLength + 0.45 * cc.MinorAxisLength)/2;
        myZeros((x.^2 + y.^2 < rad^2)) = 1;
    elseif (size(bwImage,1) < size(bwImage,2))
        myZeros((x.^2/((cc.MajorAxisLength/2)^2) + y.^2/((cc.MinorAxisLength/2)^2)) < 1) = 1;
    elseif (size(bwImage,1) > size(bwImage,2))
        myZeros((x.^2/((cc.MinorAxisLength/2)^2) + y.^2/((cc.MajorAxisLength/2)^2)) < 1) = 1;
    end
    
    image = myZeros .* image;
    image = uint8(real(ifft2(ifftshift(image))));
    
end

