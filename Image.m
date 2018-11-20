classdef Image
    %IMAGE Superclass, containing all images
    %   A list of images is represented by individual images, each of which
    %   contains several edited images of itself.
    
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
        % final image containing only DNA backbones
        bwImgThinnedRemoved;
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
        
%        Array of Points where the centers of the Nukleii are, in respect
%        to the tupper left corner
        centers;
%         linear index of centers
        indexcenters;
%         Array of Radi of the Nukleii
        radii;
        
%         Elements that have an isValid flag.
        purged;
%         array of bool values if the Nukleo I is inside the I's
%         ConnectedComponent
        contains;
        
%         x,y Dimensions of the original Image
        imgSize;
%         Array of DNA indices
        dnaNumbers
%         Array of Center of Masses of all DNA Objects, Free and Bound
        coordConnected;
%         Array of numbers of connected DNA strand for Nukleo at position i
        attachedDNA
%         array of bounding box of every DNA strand
        boundingBoxDna
%         Threshold of the first Threshhold Algo.
        thresh;
        thresh1;
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

