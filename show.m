
 
    imshowpair(imageList{index}.rawImage, imageList{index}.bwImgThickDna);
    hold on
    viscircles(imageList{index}.centers, imageList{index}.radii);
    hold on
    rad = ones(length(imageList{index}.region),1);
    viscircles(imageList{index}.region, rad);
