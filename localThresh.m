function [dnaList] = localThresh(dnaList)

for i = 1:length(dnaList)
    curr_img = dnaList{i}.image;
    local_t = threshold('otsu', curr_img);
    curr_bw = im2bw(curr_img, local_t);
    figure(1);
    subplot(1,3,1)
    imshow(curr_img);
    subplot(1,3,2)
    imshow(dnaList{i}.bwImage);
    subplot(1,3,3)
    imshow(curr_bw);
    k = waitforbuttonpress();
    
end

% function end
end