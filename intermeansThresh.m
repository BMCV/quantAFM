% This implementation is copied from OCTAVE version of graythresh
% requires input image to be of uint8
%% returns uint8 threshold
function Tout = intermeansThresh (img)
  %img = single(img);
  T = floor (mean (img (:)));
  y = hist (img(:), 0:255);
  %y = hist (img(:), 0:intmax (class (img)));

  n = numel (y) - 1;

  Tprev = NaN;

  % The threshold is found iteratively. In each iteration, the means of the
  % pixels below (mu) the threshold and above (nu) it are found. The
  % updated threshold is the mean of mu and nu.
  sumY = sum (y);
  sumB = partial_sumB (y, n);
  while T ~= Tprev
    sumAT = partial_sumA (y, T);
    sumBT = partial_sumB (y, T);

    mu = sumBT/sumAT;
    nu = (sumB-sumBT)/(sumY-sumAT);
    Tprev = T;
    T = floor((mu+nu)/2);
  end
  Tout = T;
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