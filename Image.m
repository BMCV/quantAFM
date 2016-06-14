classdef Image
    %IMAGE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        freeDnaCount int32;
        mutantDnaCount int32;
        normalDnaCount int32;
        dnaList;
        
        connectedThickDna;
        connectedThinnedDna;
        region;

        % original image as read
        rawImage;
        colorMap;
        % bw image used for artifact removal
        bwImage;      
        % first final bw image containing mostly thick DNA fragments / nucleosomes
        bwImgThickDna;
        % denoised and filtered bw image
        bwImgDen;
        % denoised, filtered and thinned bw image
        bwImgDenThinned;
        % bw image where artifacts have been filtered out
        bwFilteredImage;
        % original image after first preprocessing steps (lowpass + median filtered)
        preprocImg;
        % original image after denoising with OpenCV function
        denoisedImg;
        
        fftImage;
        fftbwImage;
        cArray;
        metaImage;
        % original image where artifacts' intensity is adjusted
        filteredImage; 
        % background of original image
        background;
        % first final bw image containing 1 px thick DNA fragments, 
        % no nucleosomes visible anymore 
        bwImgThinnedDna; 
        
        centers;
        radii;
        
        contains;
        indexcenters;
        
    end
    
    methods
        function imageObj = Image(image, colorMap)
            imageObj.freeDnaCount = 0;
            imageObj.mutantDnaCount = 0;
            imageObj.normalDnaCount = 0;
            if nargin >0
                imageObj.rawImage = image;
                imageObj.colorMap = colorMap;
            end
        end
    end
    
end

