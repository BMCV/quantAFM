function [ image ] = lowPassFilter( image )
%FFT Summary of this function goes here
%   Detailed explanation goes here

 	image = fftshift(fft2(image));
    gray_image = mat2gray(log(abs(image)));
    threshhold = graythresh(gray_image);
    bwimage = im2bw(gray_image,threshhold);
    
    SE = strel('disk', 4, 8);
    bwimage(floor(size(bwimage, 1)/2) - 20 : floor(size(bwimage, 1)/2) + 20, ...
    floor(size(bwimage, 2)/2) - 20 : floor(size(bwimage, 2)/2) + 20) = 1;
    bwimage = imerode(bwimage, SE);
    bwimage = bwareafilt(bwimage, 1, 'largest');
    bwimage = imfill(bwimage, 'holes');
    % delete the pixels in every row that has less than 75 white pixels
%     image(sum(image) < 75, :) = 0;
%    for i = 1:size(bwimage,1)
%        if(sum(bwimage(i,:)) < 75)
%            bwimage(i,:) = 0;
%        end
%    end
    % faster:
    % sum(img, 2) sums rows, returns [1,n] array
    % ~(sum(img,2) >= 75) returns inverted binary mask of where that array
    %   is larger than 75
    % find(a) returns indices of non-zero elements
    bwimage(find(~(sum(bwimage, 2)>=75)), :) = 0;
    
    SE = strel('disk', 80, 8);
    bwimage = imdilate(bwimage, SE); 
    cc = regionprops(bwimage, 'Centroid', 'MajorAxisLength', 'MinorAxisLength');
    rad = (0.55 * cc.MajorAxisLength + 0.45 * cc.MinorAxisLength)/2;
%    imshow(bwimage);
%    hold on;
%    viscircles(cc.Centroid, rad);
    xcoord = (1:size(bwimage, 1)) - cc.Centroid(1);
    ycoord = (1:size(bwimage, 2)) - cc.Centroid(2);
    
    [x,y] = ndgrid(xcoord, ycoord);
    myZeros = zeros(size(bwimage));
    myZeros((x.^2 + y.^2 < rad^2)) = 1;
    image = myZeros .* image;
    image = uint8(real(ifft2(ifftshift(image))));
    
end

