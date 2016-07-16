function [ bool ] = isLPoint( picture )
% picture is a 3x3 matrix with the potential LPoint in the middle
% returns if this point is a LPoint or not

% starting assumption is that it is no LPoint
bool = 0;

    if isequal(picture, [0,1,0;1,1,0;0,0,0])
        bool = 1;
    end
    
    if isequal(picture, [0,0,0;1,1,0;0,1,0])
        bool = 1;
    end
    
    if isequal(picture, [0,1,0;0,1,1;0,0,0])
        bool = 1;
    end
    
    if isequal(picture, [0,0,0;0,1,1;0,1,0])
        bool = 1;
    end
end

