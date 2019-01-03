function showFT(im_ft)
    imshow(mat2gray(100*log(1+abs(fftshift(im_ft)))),[])