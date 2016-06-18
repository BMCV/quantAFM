figure;
imshow(imageList{index}.dnaList{dnaIndex}.bwImageThinned);
if(~strcmp(imageList{index}.dnaList{dnaIndex}.type , 'free'))
text(20,1,num2str(imageList{index}.dnaList{dnaIndex}.angle1));
text(20,20,num2str(imageList{index}.dnaList{dnaIndex}.angle2));
end