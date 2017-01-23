function done = writeToCsvFile(filename, imageObj)

    global PIXELLENGTH REALVALUE;
    dnaCount = size(imageObj.dnaList,2);
    number = zeros(dnaCount,1) -1;
    xm = zeros(dnaCount,1) -1;
    ym = zeros(dnaCount,1) -1;
    radius = zeros(dnaCount,1) -1;
    length = zeros(dnaCount,1) -1;
    length_botharms = zeros(dnaCount,1) -1;
    short_arm = zeros(dnaCount,1) -1;
    long_arm = zeros(dnaCount,1) -1;
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
        if (curr.hasNucleus == 0 ) % no nucleosome detected
            angle1(dnaIndex,1) = 0;
            angle2(dnaIndex,1) = 0;
            radius(dnaIndex,1) = 0;
            length(dnaIndex,1) = curr.length{1};
            length_botharms(dnaIndex,1) = 0;
            short_arm(dnaIndex,1) = 0;
            long_arm(dnaIndex,1) = 0;
            numNucleosomes(dnaIndex,1) = 0;
        else  % at least one nucleosome present
           angle1(dnaIndex,1) = curr.angle1;
           angle2(dnaIndex,1) = curr.angle2;
           radius(dnaIndex,1) = curr.attachedNukleo{1}.rad;
           numNucleosomes(dnaIndex,1) =  numel(curr.attachedNukleo) ;
           if size(curr.length, 2) == 3 % one short and one long arm found
                if curr.length{2} > curr.length{3}
                    short_arm(dnaIndex,1) = curr.length{3};
                    long_arm(dnaIndex,1) = curr.length{2};
                else
                    short_arm(dnaIndex,1) = curr.length{2};
                    long_arm(dnaIndex,1) = curr.length{3};
                end
                % now differ between total length and added length.
                length_botharms(dnaIndex,1) = curr.length{2} + curr.length{3};
%                 length(dnaIndex,1) = curr.length{2}+curr.length{3};
                length(dnaIndex,1) = curr.length{1};
           else
               short_arm(dnaIndex,1) = 0;
               long_arm(dnaIndex,1) = curr.length{1};
               %shouldn't long arm be 0?
%                long_arm(dnaIndex,1) = 0;
               length(dnaIndex,1) = curr.length{1};
           end 
        end
    end
    
    if (REALVALUE == 1)
       length = length * PIXELLENGTH;
       length_botharms = length_botharms * PIXELLENGTH;
       short_arm = short_arm * PIXELLENGTH;
       long_arm = long_arm * PIXELLENGTH;
       radius = radius * PIXELLENGTH;
    end
    
    T = table(number, xm, ym, length, hasNucleus, length_botharms, short_arm, long_arm ,...
        radius, isValid, angle1, angle2, numNucleosomes);
    writetable(T, filename);
    done = 'done';
end