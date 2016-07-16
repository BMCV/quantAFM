figure;
imshow(imageList{index}.dnaList{dnaIndex}.bwImageThinned);


text(imageList{index}.dnaList{dnaIndex}.position(:,1),...
    imageList{index}.dnaList{dnaIndex}.position(:,2), ...
    num2str(imageList{index}.dnaList{dnaIndex}.number), 'Color' , 'red');
