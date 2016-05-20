function  t  = threshold( method , image )
%T Summary of this function goes here
%   Detailed explanation goes here
          
    switch method
        case 'otsu'
            t = graythresh(image);
        case 'intermeans'
            t = intermeansThresh(image)/255;
        case 'moments'
            t = momentsThresh(image)/255;
        otherwise 
            t = 0;
            disp('Unknown Threshold parameter "name"');
    end


            
         

end

