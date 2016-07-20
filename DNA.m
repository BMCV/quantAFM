classdef (Abstract) DNA
    %DNA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
%         Upper boundary for DNA length
        MAX_LENGTH = 200;
%         Lower boundary for DNA length
        MIN_LENGTH = 20;
%  Number of DNA strand
        number;
%  Length of Dna strand
        length;
%        flag that identifies whether DNA object is valid or not;
%        a 1 represents a valid object, a 0 a non-valid object
%        DNA object is invalid if it is
%           - too short (< 20 pixel)
%           - too long (> 200 pixel)
%           - ... (more criteria?)
        isValid;
%         center of mass of DNA strand, respect to upper left corner
        position = zeros(1,2) ; %center of mass
%        Boolean, is 1 if DNA has nucleus, else 0
        hasNucleus;
%       Array of linear Indices of the DNA strand in the bwDnaThickImg 
        connectedThick;
%        Array of linear indices of the DNA strand in the bwDnaThinImg
        connectedThinned;
%        Array of linear indices of the DNA backbone in the bwDnaThinImg
        connectedThinnedRemoved;
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


