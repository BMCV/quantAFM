
 
%     imshowpair(imageList{index}.rawImage, imageList{index}.bwImgThickDna);
    imshow(label2rgb(labelmatrix(imageList{index}.connectedThinnedDna)))
    hold on
    viscircles(imageList{index}.centers, imageList{index}.radii);
    hold on
    rad = ones(length(imageList{index}.region),1);
    viscircles(imageList{index}.region, rad);
%     text(imageList{index}.region(:,1), imageList{index}.region(:,2),...
%         [repmat( ' ', size(imageList{index}.region,1)),num2str((1:1:size(imageList{index}.region,1))')],'Color', 'b') ;
% text(imageList{index}.region(:,1), imageList{index}.region(:,2),num2str([1:size(imageList{index}.region,1)]'));
%text(imageList{index}.region(:,1), imageList{index}.region(:,2), num2str ([1: size(imageList{index}.region,1)]'));

% text(imageList{index}.centers(:,1), imageList{index}.centers(:,2), num2str ([1: size(imageList{index}.centers,1)]'));
text(imageList{index}.centers(:,1), imageList{index}.centers(:,2), num2str(imageList{index}.attachedDNA));
