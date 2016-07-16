img = zeros(dna.sizeImg);
img(dna.bwImage)=1.;
img(mask) = 0.75;
img(dnaImgThinned) = 0.5;
figure; imshow(img);
hold on;
viscircles(dna.attachedNukleo{1}.localCenter, dna.attachedNukleo{1}.rad);
%scatter(coords(:,1),coords(:,2),'.', 'green')
scatter(p(:,1),p(:,2), 'x','blue')
plot([p(1,1) dna.attachedNukleo{1}.localCenter(1)],[p(1,2) dna.attachedNukleo{1}.localCenter(2)],'blue')
plot([p(2,1) dna.attachedNukleo{1}.localCenter(1)],[p(2,2) dna.attachedNukleo{1}.localCenter(2)],'blue')
