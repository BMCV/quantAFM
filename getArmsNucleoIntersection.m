function arms = getArmsNucleoIntersection(dnaObj, short)
    arms = 0;
    % cut out nucleo
    if short == 0
        dnaImgThinned = dnaObj.bwImageThinnedRemoved;
    else
        dnaImgThinned = dnaObj.bwImageThinnedRemovedShort;
    end
    mask = bsxfun(@plus, ((1:dnaObj.sizeImg(2)) - ...
        dnaObj.attachedNukleo{1}.localCenter(1)).^2, ...
        ((transpose(1:dnaObj.sizeImg(1)) - ...
        dnaObj.attachedNukleo{1}.localCenter(2)).^2)) ...
            < dnaObj.attachedNukleo{1}.rad^2;
    dnaImgThinned(mask) = 0;
    % Get the two branches
    arms = bwconncomp(dnaImgThinned);
end