function done = writeToCsvFile(filename, imageObj)
    dnaCount = size(imageObj.dnaList,2);
    number = {};
    xm = {};
    ym = {};
    length = {};
    type = {};
    angle1 = {};
    angle2 = {};
    numNucleosomes = {};
    for dnaIndex = 1:dnaCount
        curr = imageObj.dnaList{dnaIndex};
        number{dnaIndex} = curr.number;
        xm{dnaIndex} = curr.position(1);
        ym{dnaIndex} = curr.position(2);
        length{dnaIndex} = curr.length;
        type{dnaIndex} = curr.type;
        if (strcmp(curr.type,'free'))
            angle1{dnaIndex} = 0;
            angle2{dnaIndex} = 0;
            numNucleosomes{dnaIndex} = 0;
        else
           angle1{dnaIndex} = curr.angle1;
           angle2{dnaIndex} = curr.angle2;
           numNucleosomes{dnaIndex} = size(curr.attachedNukleo, 2); 
        end 
    end
    number = number';
    xm = xm'; 
    ym = ym';
    length = length';
    type = type';
    angle1 = angle1';
    angle2 = angle2';
    numNucleosomes = numNucleosomes';
    T = table(number, xm, ym, length, type, angle1, angle2, numNucleosomes);
    writetable(T, filename);
    done = 'done';
end