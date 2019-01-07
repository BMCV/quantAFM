function done = showAngle3(dna, mask, mask2, radius)
    figure; 
    subplot(1,2,1);
    img = zeros(dna.sizeImg);
    img(dna.bwImage)=1.;
    img(mask) = 0.75;
    img(dna.bwImageThinned) = 0.5;
    imshow(img);
    hold on;
    viscircles( dna.attachedNukleo{1}.localCenter, dna.attachedNukleo{1}.rad);
    
    subplot(1,2,2);
    img = zeros(dna.sizeImg);
    img(dna.bwImage)=1.;
    img(mask2) = 0.75;
    img(dna.bwImageThinned) = 0.5;
    imshow(img);
    hold on;
    viscircles( dna.attachedNukleo{1}.localCenter, radius);
    w = waitforbuttonpress();
end