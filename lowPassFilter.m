function [ image ] = lowPassFilter( image )
%this function gets an non bw image and uses Fast Fourier Transformation to
%filter high frequencies. The output is the filtered image.
    
    % transform image to Fourier space and transforms the magnitude image
    % into a binary image
 	image = fftshift(fft2(image));
    gray_image = mat2gray(log(abs(image)));
    threshold = graythresh(gray_image);
    bwImage = im2bw(gray_image,threshold);
    
    % In the magnitude image one can see that we have a brighter area in
    % the middle of the picture. We now want to fit a circle into the
    % picture that approximates the white area in the middle and separates
    % it from the lower frequencies outside of the circle
    SE = strel('disk', 4, 8);
    % color a little rectangle in the middle of the image white to avoid
    % some edge cases that create errors when there is too much black in
    % the middle
    bwImage(floor(size(bwImage, 1)/2) - 20 : floor(size(bwImage, 1)/2) + 20, ...
    floor(size(bwImage, 2)/2) - 20 : floor(size(bwImage, 2)/2) + 20) = 1;
    % erode image so that all white pixels outside of the center get
    % deleted and only the space in the middle of the image remains
    bwImage = imerode(bwImage, SE);
    % take only the largest connected area that remains
    bwImage = bwareafilt(bwImage, 1, 'largest');
    % fill holes in this area
    bwImage = imfill(bwImage, 'holes');

    % this deletes all white pixels in rows that have less than 75 pixels
    % so the white line from top to bottom gets removed and only the circle
    % in the middle remains
    bwImage(find(~(sum(bwImage, 2)>=75)), :) = 0;
    
    % now the little circle in the middle of the image gets dilated since
    % we need to revert the eroding from earlier. we do this in a way so
    % the area we receive is representative for our image
    SE = strel('disk', 80, 8);
    bwImage = imdilate(bwImage, SE);
    % get properties we need to approximate our circle
    cc = regionprops(bwImage, 'Centroid', 'MajorAxisLength', 'MinorAxisLength');

    % create grid so that in the picture we create the center of our circle
    % has the value 0 and increases linearly to the edges
    xcoord = (1:size(bwImage, 2)) - cc.Centroid(1);
    ycoord = (1:size(bwImage, 1)) - cc.Centroid(2);
    
    [y,x] = ndgrid(ycoord, xcoord);
    % create image with size of our original image only filled with 0s
    myZeros = zeros(size(bwImage));
    
    % now fit a circle (or ellipse for non quadratic pictures) into the
    % myZeros image. for that we use the grid we created which has lower
    % values near the centre of the circle (and gets linearly higher the
    % farer we get from the center).
    if (size(bwImage,1) == size(bwImage,2))
        % for circles we use the major and minor axis lentgh of our area
        rad = (0.55 * cc.MajorAxisLength + 0.45 * cc.MinorAxisLength)/2;
        myZeros((x.^2 + y.^2 < rad^2)) = 1;
    elseif (size(bwImage,1) < size(bwImage,2))
        % ellipse equation
        myZeros((x.^2/((cc.MajorAxisLength/2)^2) + y.^2/((cc.MinorAxisLength/2)^2)) < 1) = 1;
    elseif (size(bwImage,1) > size(bwImage,2))
        myZeros((x.^2/((cc.MinorAxisLength/2)^2) + y.^2/((cc.MajorAxisLength/2)^2)) < 1) = 1;
    end
    
    % in my Zeros is now our fitted circle. all other values are 0. so we
    % now delete all values in image that are outside of the circle by
    % multiplying them.
    image = myZeros .* image;
    % backtransformation from frequency space
    image = uint8(real(ifft2(ifftshift(image))));
    
end

