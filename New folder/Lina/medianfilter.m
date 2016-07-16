function im_f = medianfilter(im)
    %im_f = medfilt2(im);
    % http://stackoverflow.com/questions/27535535/matlab-median-filter-code
    im_pad = padarray(im, [1 1]);
    im_col = im2col(im_pad, [3 3], 'sliding');
    sorted_cols = sort(im_col, 1, 'ascend');
    im_f = col2im(sorted_cols(5, :), [3 3], size(im_pad), 'sliding');
    im_f = mat2gray(im_f);
    figure('Name', 'Median filter')
    subplot(2,2,1)
    imshow(im)
    title('Original image')
    subplot(2,2,2)
    imshow(im_f)
    title('Median filtered image')

    % Fourier Analysis
    ft_im =fft2(im);
    %ft_im = mat2gray(ft_im);
    subplot(2,2,3)
    showFT(ft_im)
    title('FT original image')
    % Fourier Analysis
    ft_bg =fft2(im_f);
    %ft_bg = mat2gray(ft_bg);
    subplot(2,2,4)
    showFT(ft_bg)
    title('FT Median filtered image')