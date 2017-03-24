function done = export_pixels(filename, imageObj)
%EXPORT_PIXELS Outputs the position of pixels to csv/tsv files
%   To calculate the stiffness of the DNA molecules, specific data formats
%   are required by the appropriate program.
%   To ensure the compatibility, the user may enable or disable to
%   conversion to the absolute image position (depending on the scansize of
%   the original image).

    global PIXELLENGTH REALVALUE;
    dnaCount = size(imageObj.dnaList,2);
    
    for i = 1:dnaCount
        
        
    end

end

