% @Author:  Linbo<linbo.me>
% @Version: 1.0  25/10/2014
% This is the implementation of the 
% Zhang-Suen Thinning Algorithm for skeletonization.

function BW_Thinned = thinningZhangSuen(BW_Original)
pad = 1;
changing = 1;
BW_Thinned = padarray(BW_Original,[pad,pad]);
[rows, columns] = size(BW_Thinned);
BW_Del = ones(rows, columns);
while changing
    % BW_Del = ones(rows, columns);
    changing = 0;
    % Step 1
    for i=2:rows-1
        for j = 2:columns-1
            P = [BW_Thinned(i,j) BW_Thinned(i-1,j) BW_Thinned(i-1,j+1) BW_Thinned(i,j+1) BW_Thinned(i+1,j+1) BW_Thinned(i+1,j) BW_Thinned(i+1,j-1) BW_Thinned(i,j-1) BW_Thinned(i-1,j-1) BW_Thinned(i-1,j)]; % P1, P2, P3, ... , P8, P9, P2
            if (BW_Thinned(i,j) == 1 &&  sum(P(2:end-1))<=6 && sum(P(2:end-1)) >=2 && P(2)*P(4)*P(6)==0 && P(4)*P(6)*P(8)==0)   % conditions
                % No. of 0,1 patterns (transitions from 0 to 1) in the ordered sequence
                A = 0;
                for k = 2:size(P,2)-1
                    if P(k) == 0 && P(k+1)==1
                        A = A+1;
                    end
                end
                if (A==1)
                    BW_Del(i,j)=0;
                    changing = 1;
                end
            end
        end
    end
    BW_Thinned = BW_Thinned.*BW_Del;  % the deletion after all the pixels have been visited
    % Step 2 
    for i=2:rows-1
        for j = 2:columns-1
            P = [BW_Thinned(i,j) BW_Thinned(i-1,j) BW_Thinned(i-1,j+1) BW_Thinned(i,j+1) BW_Thinned(i+1,j+1) BW_Thinned(i+1,j) BW_Thinned(i+1,j-1) BW_Thinned(i,j-1) BW_Thinned(i-1,j-1) BW_Thinned(i-1,j)];
            if (BW_Thinned(i,j) == 1 && sum(P(2:end-1))<=6 && sum(P(2:end-1)) >=2 && P(2)*P(4)*P(8)==0 && P(2)*P(6)*P(8)==0)   % conditions
                A = 0;
                for k = 2:size(P,2)-1
                    if P(k) == 0 && P(k+1)==1
                        A = A+1;
                    end
                end
                if (A==1)
                    BW_Del(i,j)=0;
                    changing = 1;
                end
            end
        end
    end
    BW_Thinned = BW_Thinned.*BW_Del;
end%while
BW_Thinned = BW_Thinned(1+pad:end-pad,1+pad:end-pad);
%figure
%subplot(1,2,1)
%imshow(BW_Original)
%subplot(1,2,2)
%imshow(BW_Thinned)
