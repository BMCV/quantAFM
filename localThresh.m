function [dnaList, average, average_skeleton] = localThresh(dnaList)
% currently not in use. Test for a local (second) thresholding mechanism
% that prevents objects with small angles (i.e. probable
% self-intersections) to touch in the binarized image.
average = zeros(length(dnaList),1);
average_skeleton = zeros(length(dnaList),1);

for i = 1:length(dnaList)
    % create local thresholded image
    curr_img = dnaList{i}.image;
    local_t = threshold('otsu', curr_img);
    curr_bw = im2bw(curr_img, local_t);
    
    % calculate averages
    temp_gray = uint8(dnaList{i}.bwImage) .* dnaList{i}.image;
    temp_local_gray = uint8(curr_bw) .* dnaList{i}.image;
    temp_skeleton = uint8(dnaList{i}.bwImageThinned) .* dnaList{i}.image;
    average(i) = sum(sum(temp_gray)) / sum(sum(dnaList{i}.bwImage));
    average_skeleton(i) = sum(sum(temp_skeleton)) / sum(sum(dnaList{i}.bwImageThinned));
    
    local_average = sum(sum(temp_local_gray)) / sum(sum(curr_bw));
    % create comparison image with removed skeltal values.
    removed_skeleton = logical(temp_skeleton > average(i));
    removed_local_skeleton = logical(temp_skeleton > local_average);
    
    figure(1);
    subplot(1,6,1)
    imshow(curr_img);
    subplot(1,6,2)
    imshow(dnaList{i}.bwImage);
    subplot(1,6,3)
    imshow(curr_bw);
    subplot(1,6,4)
    imshow(dnaList{i}.bwImageThinned);
    subplot(1,6,5)
    imshow(removed_skeleton);
    subplot(1,6,6)
    imshow(removed_local_skeleton);
    k = waitforbuttonpress();
    
end

% function end
end