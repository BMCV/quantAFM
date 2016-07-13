
 
%     imshowpair(imageList{index}.rawImage, imageList{index}.bwImgThickDna);
    %imshow(label2rgb(labelmatrix(imageList{index}.connectedThinnedDna));
figure;
imshow(imfuse(imageList{index}.rawImage, imageList{index}.bwImgThickDna));
hold on
viscircles(imageList{index}.centers, imageList{index}.radii, ...
    'LineWidth', 0.5, 'EnhanceVisibility', false);
hold on
rad = ones(length(imageList{index}.region),1);
viscircles(imageList{index}.region, rad,...
    'LineWidth', 0.5, 'EnhanceVisibility', false);
%     text(imageList{index}.region(:,1), imageList{index}.region(:,2),...
%         [repmat( ' ', size(imageList{index}.region,1)),num2str((1:1:size(imageList{index}.region,1))')],'Color', 'b') ;
% text(imageList{index}.region(:,1), imageList{index}.region(:,2),num2str([1:size(imageList{index}.region,1)]'));
%text(imageList{index}.region(:,1), imageList{index}.region(:,2), num2str ([1: size(imageList{index}.region,1)]'));

 text(imageList{index}.region(:,1), imageList{index}.region(:,2), ...
     num2str ([1: size(imageList{index}.region,1)]'), 'FontSize', 6);
 %text(imageList{index}.centers(:,1), imageList{index}.centers(:,2), num2str(imageList{index}.attachedDNA));
 hold off;
