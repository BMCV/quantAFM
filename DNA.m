classdef (Abstract) DNA
    %DNA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
%  Number of DNA strand
        number;
%  Length of Dna strand
        length ;
        
%         center of mass of DNA strand, respect to upper left corner
        position = zeros(1,2) ; %center of mass
%        String Name of type
        type char;
%       Array of linear Indices of the DNA strand in the bwDnaThickImg 
        connectedThick;
%        Array of linear indices of the DNA strand in the bwDnaThinImg
        connectedThinned;
%         BW Image
        bwImage;
%         BWimage Thin DNA
        bwImageThinned;
%         bwImage with removes arms
        bwImageThinnedRemoved
%         size of the subImage/bwImage
        sizeImg
%         Bounding Box struct
        bBox;
    end
    
    methods
    end
    
end


