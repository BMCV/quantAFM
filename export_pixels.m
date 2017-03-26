function done = export_pixels(folder, imageObj)
%EXPORT_PIXELS Outputs the position of pixels to csv/tsv files
%   To calculate the stiffness of the DNA molecules, specific data formats
%   are required by the appropriate program.
%   To ensure the compatibility, the user may enable or disable to
%   conversion to the absolute image position (depending on the scansize of
%   the original image).

    global PIXELLENGTH EXPORTREAL EXPORTONLYVALID;
    dnaCount = size(imageObj.dnaList,2);
    
    for i = 1:dnaCount
        if ((EXPORTONLYVALID == 1 && imageObj.dnaList{i}.isValid ==1) || (EXPORTONLYVALID == 0))
            [x,y] = ind2sub(imageObj.dnaList{i}.sizeImg, imageObj.dnaList{i}.connectedThinnedRemoved);
            if (EXPORTREAL==1)
                x = round(x*PIXELLENGTH, 6);
                y = round(y*PIXELLENGTH, 6);
            end
            output = table(x,y);
            
            temp = [folder 'file' num2str(i, '%03i') '.txt'];
            writetable(output, temp, 'delimiter', ' ', 'WriteVariableNames', 0)
        end
    end

    done = 'done';
end

