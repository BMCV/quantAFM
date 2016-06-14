% This implementation is copied from OCTAVE version of graythresh
% requires input image to be of uint8
%% returns uint8 threshold
function level = momentsThresh (img)
  y = hist (img(:), 0:255);
  n = numel (y) - 1;

  %% The threshold is chosen such that partial_sumA(y,t)/partial_sumA(y,n)
  %% is closest to x0.
  sumY = sum (y);
  Avec = (cumsum (y)- y(1)) / sumY;

  sumB = partial_sumB (y,n);
  sumC = partial_sumC (y,n);
  sumD = partial_sumD (y,n);
  %% The following finds x0.
  x2 = (sumB*sumC - sumY*sumD) / (sumY*sumC - sumB^2);
  x1 = (sumB*sumD - sumC^2) / (sumY*sumC - sumB^2);
  x0 = .5 - (sumB/sumY + x2/2) / sqrt (x2^2 - 4*x1);

  %% And finally the threshold
  [~, ind] = min (abs (Avec-x0));
  level = ind-1;
end


%% partial sums from C. A. Glasbey, "An analysis of histogram-based thresholding
%% algorithms," CVGIP: Graphical Models and Image Processing, vol. 55, pp. 532-537, 1993.
function x = partial_sumA (y, j)
  x = sum (y(1:j+1));
end
function x = partial_sumB (y, j)
  ind = 0:j;
  x   = ind*y(1:j+1)';
end
function x = partial_sumC (y, j)
  ind = 0:j;
  x = ind.^2*y(1:j+1)';
end
function x = partial_sumD (y, j)
  ind = 0:j;
  x = ind.^3*y(1:j+1)';
end
