figure;
imshow(imageList{index}.dnaList{dnaIndex}.bwImageThinned);
if(~strcmp(imageList{index}.dnaList{dnaIndex}.type , 'free'))
text(1,1,num2str(imageList{index}.dnaList{dnaIndex}.angle1),'Color' , 'red');
text(20,1,num2str(imageList{index}.dnaList{dnaIndex}.angle2),'Color' , 'red');
text(imageList{index}.dnaList{dnaIndex}.attachedNukleo{1}.localCenter(:,1),...
    imageList{index}.dnaList{dnaIndex}.attachedNukleo{1}.localCenter(:,2), ...
    num2str(imageList{index}.dnaList{dnaIndex}.number), 'Color' , 'red');
end