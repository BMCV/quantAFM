function im_bin = thresh_double(im_bin,thresh)
    figure
    subplot(2,2,1)
    imshow(im_bin)
    title('Input Image')
    im_bin(im_bin < thresh) = 0;
    subplot(2,2,2)
    imshow(im_bin)
    title('Simple Threshold')
    
    im_bin = bwmorph(im_bin,'thin',Inf);
    subplot(2,2,3)
    imshow(im_bin)
    title('Thinned')
    [labeledImage, numRegions] = bwlabel(im_bin);
    for i=1:numRegions
        [r, c] = find(labeledImage==i);
        if (size(r,1)<20 || nnz(r==1|r==size(im_bin,1)|c==1|c==size(im_bin,2)) )
            im_bin(r,c) = 0;
        end
    end

    subplot(2,2,4)
    imshow(im_bin)
    title('Removed Noise')