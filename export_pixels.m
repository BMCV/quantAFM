function done = export_pixels(filename, imageObj)
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
                x = x*PIXELLENGTH;
                y = y*PIXELLENGTH;
            end
            output = table(x,y);
            
            temp = filename
            if( strcmp(getenv('OS'),'Windows_NT'))
                
                temp = [temp '_' num2str(i, '%02i') '.csv'];
            else
                temp = [temp '_' num2str(i, '%02i') '.csv'];
            end
            writetable(output, temp, 'delimiter', ' ', 'WriteVariableNames', 0)
        end
    end

end

