classdef Image
    %IMAGE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        freeDnaCount int32;
        mutantDnaCount int32;
        normalDnaCount int32;
        DnaMap = containers.Map
        
        connected;
        region;
        
        image;
        colorMap;
        bwImage;
        bwImageremoved;
        medImage;
        fftImage;
        fftbwImage;
        cArray;
        metaImage;
        complementImage;
        
        
    end
    
    methods
        function imageObj = Image(image, colorMap)
            imageObj.freeDnaCount = 0;
            imageObj.mutantDnaCount = 0;
            imageObj.normalDnaCount = 0;
            if nargin >0
                imageObj.image = image;
                imageObj.colorMap = colorMap;
            end
        end
        function addMapPair(dna)
            key = strcat(num2hex(dna.position(1)), num2hex(dna.positon(2)));
            this.DnaMap(key) = dna;
        end
    end
    
end

