function done = writeToCsvFile(filename, imageObj)
    dnaCount = size(imageObj.dnaList,2);
    number = zeros(dnaCount,1) -1;
    xm = zeros(dnaCount,1) -1;
    ym = zeros(dnaCount,1) -1;
    radius = zeros(dnaCount,1) -1;
    length = zeros(dnaCount,1) -1;
    length_arm1 = zeros(dnaCount,1) -1;
    length_arm2 = zeros(dnaCount,1) -1;
    hasNucleus = zeros(dnaCount,1) -1;
    isValid = zeros(dnaCount,1) -1;
    angle1 = zeros(dnaCount,1) -1;
    angle2 = zeros(dnaCount,1) -1;
    numNucleosomes = zeros(dnaCount,1) -1;
    for dnaIndex = 1:dnaCount
        curr = imageObj.dnaList{dnaIndex};
        number(dnaIndex,1) = curr.number;
        xm(dnaIndex,1) = curr.position(2);
        ym(dnaIndex,1) = curr.position(1);
        hasNucleus(dnaIndex,1) = curr.hasNucleus;
        isValid(dnaIndex,1) = curr.isValid;
        if (curr.hasNucleus == 0 || numel(curr.attachedNukleo) ~= 1)
            angle1(dnaIndex,1) = 0;
            angle2(dnaIndex,1) = 0;
            radius(dnaIndex,1) = 0;
            length(dnaIndex,1) = curr.length{1};
            length_arm1(dnaIndex,1) = 0;
            length_arm2(dnaIndex,1) = 0;
            numNucleosomes(dnaIndex,1) = 0;
        else
           angle1(dnaIndex,1) = curr.angle1;
           angle2(dnaIndex,1) = curr.angle2;
           radius(dnaIndex,1) = curr.attachedNukleo{1}.rad;
           numNucleosomes(dnaIndex,1) = size(curr.attachedNukleo, 2);
           length(dnaIndex,1) = 0;
           length_arm1(dnaIndex,1) = curr.length{1};
           if size(curr.length, 2) == 2
               length_arm2(dnaIndex,1) = curr.length{2};
           else
               length_arm2(dnaIndex,1) = 0;
           end 
    end
    T = table(number, xm, ym, length, hasNucleus, length_arm1, length_arm2 ,...
        radius, isValid, angle1, angle2, numNucleosomes);
    writetable(T, filename);
    done = 'done';
end