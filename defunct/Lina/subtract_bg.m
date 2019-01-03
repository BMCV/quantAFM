function im_bg = subtract_bg(im);
    kernel = [-1 -2 -1; -2 12 -2; -1 -2 -1]/16;
    filteredImage = imfilter(single(im), kernel,'replicate','same');
    filteredImage = mat2gray(filteredImage);
    im_bg = mat2gray(im-filteredImage);
    figure%('Name', 'Background_Subtraction')
    subplot(1,3,1)
    imshow(im)
    title('Original Image')
    subplot(1,3,2)
    imshow(filteredImage, []);
    title('Filtered Image');
    subplot(1,3,3)
    imshow(im_bg);
    title('Subtracted BG Image')