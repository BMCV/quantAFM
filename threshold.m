function  t  = threshold( method , img )
%T Summary of this function goes here
%   Detailed explanation goes here
    ihist = hist (img(:), 0:255) ;
    
    switch method
        case 'otsu'
            t = graythresh(img); % Matlab implementation
        case 'intermeans'
            t = intermeansThresh(img)/255;
        case 'moments'
            t = momentsThresh(ihist)/255;
        case 'concavity'
            t = concavity(img)/255;
        case 'intermodes'
            t = intermodes(ihist)/255;
        case 'maxentropy'
            t = maxentropy(ihist)/255;
        case 'maxlikelihood'
            t = maxlikelihood(ihist)/255;
        case 'minerror'
            t = minerror_iter(ihist)/255;
        case 'minimum'
            t = minimum(img)/255;
        case 'percentile'
            t = percentile(ihist)/255;
        otherwise 
            t = 0;
            disp('Unknown Threshold parameter "name"');
    end

end

%% All following thesholding functions are modified from Octave's 
%% graythesh method. Code to be found at:
%% https://sourceforge.net/p/octave/image/ci/default/tree/inst/graythresh.m#l175
% Copyright (C) 2004 Antti Niemistö <antti.niemisto@tut.fi>
% Copyright (C) 2007 Søren Hauberg
% Copyright (C) 2012-2013 Carnë Draug <carandraug+dev@gmail.com>
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, see <http://www.gnu.org/licenses/>.
% @item Otsu (uses Matlab function)
% Implements Otsu's method as described in @cite{Nobuyuki Otsu (1979). "A
% threshold selection method from gray-level histograms", IEEE Trans. Sys.,
% Man., Cyber. 9 (1): 62-66}. This algorithm chooses the threshold to minimize
% the intraclass variance of the black and white pixels.
%
% The second output, @var{sep} represents the ``goodness'' (or separability) of
% the threshold at @var{level}.  It is a value within the range [0 1], the
% lower bound (zero) being attainable by, and only by, histograms having a
% single constant gray level, and the upper bound being attainable by, and only
% by, two-valued pictures.
%
% @item concavity
% Find a global threshold for a grayscale image by choosing the threshold to
% be in the shoulder of the histogram @cite{A. Rosenfeld, and P. De La Torre
% (1983). "Histogram concavity analysis as an aid in threshold selection", IEEE
% Transactions on Systems, Man, and Cybernetics, 13: 231-235}.
%
% @item intermodes
% This assumes a bimodal histogram and chooses the threshold to be the mean of
% the two peaks of the bimodal histogram @cite{J. M. S. Prewitt, and M. L.
% Mendelsohn (1966). "The analysis of cell images", Annals of the New York
% Academy of Sciences, 128: 1035-1053}.
%
% Images with histograms having extremely unequal peaks or a broad and flat
% valley are unsuitable for this method.
%
% @item intermeans
% Iterative procedure based on the iterative intermeans algorithm of @cite{T.
% Ridler, and S. Calvard (1978). "Picture thresholding using an iterative
% selection method", IEEE Transactions on Systems, Man, and Cybernetics, 8: 630-632}
% and @cite{H. J. Trussell (1979). "Comments on 'Picture thresholding using an
% iterative selection method'", IEEE Transactions on Systems, Man, and Cybernetics,
% 9: 311}.
%
% Note that several implementations of this method exist. See the source code
% for details.
%
% @item MaxEntropy
% Implements Kapur-Sahoo-Wong (Maximum Entropy) thresholding method based on the
% entropy of the image histogram @cite{J. N. Kapur, P. K. Sahoo, and A. C. K. Wong
% (1985). "A new method for gray-level picture thresholding using the entropy
% of the histogram", Graphical Models and Image Processing, 29(3): 273-285}.
%
% @item MaxLikelihood
% Find a global threshold for a grayscale image using the maximum likelihood
% via expectation maximization method @cite{A. P. Dempster, N. M. Laird, and D. B.
% Rubin (1977). "Maximum likelihood from incomplete data via the EM algorithm",
% Journal of the Royal Statistical Society, Series B, 39:1-38}.
%
% @item mean
% The mean intensity value. It is mostly used by other methods as a first guess
% threshold.
%
% @item MinError
% An iterative implementation of Kittler and Illingworth's Minimum Error
% thresholding @cite{J. Kittler, and J. Illingworth (1986). "Minimum error
% thresholding", Pattern recognition, 19: 41-47}.
%
% This implementation seems to converge more often than the original.
% Nevertheless, sometimes the algorithm does not converge to a solution. In
% that case a warning is displayed and defaults to the initial estimate of the
% mean method.
%
% @item minimum
% This assumes a bimodal histogram and chooses the threshold to be in the
% valley of the bimodal histogram.  This method is also known as the mode
% method @cite{J. M. S. Prewitt, and M. L. Mendelsohn (1966). "The analysis of
% cell images", Annals of the New York Academy of Sciences, 128: 1035-1053}.
%
% Images with histograms having extremely unequal peaks or a broad and flat
% valley are unsuitable for this method.
%
% @item moments
% Find a global threshold for a grayscale image using moment preserving
% thresholding method @cite{W. Tsai (1985). "Moment-preserving thresholding:
% a new approach", Computer Vision, Graphics, and Image Processing, 29: 377-393}
%
% @item percentile
% Assumes a specific fraction of pixels (set at @var{options}) to be background.
% If no value is given, assumes 0.5 (equal distribution of background and foreground)
% @cite{W Doyle (1962). "Operation useful for similarity-invariant pattern
% recognition", Journal of the Association for Computing Machinery 9: 259-267}
% @end table
%
% @seealso{im2bw}
% @end deftypefn

% Notes:
%  * The following methods were adapted from http://www.cs.tut.fi/~ant/histthresh/
%      intermodes    percentile      minimum
%      MaxEntropy    MaxLikelihood   intermeans
%      moments       minerror        concavity
%  * Carnë Draug implemented and vectorized the Otsu's method
%  * Carnë Draug vectorized percentile and moments.
%  * missing methods from ImageJ
%      Yen           triangle        RenyiEntropy
%      Shanbhag      Li              Huang
%      ImageJ




%% Concavity thresholding
function T = concavity (h)
  n = numel (h) - 1;
  H = hconvhull(h);

  % Find the local maxima of the difference H-h.
  lmax = flocmax(H-h);

  % Find the histogram balance around each index.
  for k = 0:n
    E(k+1) = hbalance(h,k);
  end

  % The threshold is the local maximum with highest balance.
  E = E.*lmax;
  [~ ,ind] = max(E);
  T = ind-1;
end

%% Intermeans thresholding algorithm
% requires input image to be of uint8
% returns uint8 threshold
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

%% Intermodes thresholding
function [T] = intermodes (y)
  % checked with ImageJ and is slightly different but not by much
  n = numel (y) - 1;

  % Smooth the histogram by iterative three point mean filtering.
  iter = 0;
  while ~bimodtest(y)
    h = ones(1,3)/3;
    y = conv2(y,h,'same');
    iter = iter+1;
    % If the histogram turns out not to be bimodal, set T to zero.
    if iter > 10000;
      T{1} = 0;
      return
    end
  end

  % The threshold is the mean of the two peaks of the histogram.
  ind = 0;
  for k = 2:n
    if y(k-1) < y(k) && y(k+1) < y(k)
      ind = ind+1;
      TT(ind) = k-1;
    end
  end
  T{1} = floor(mean(TT));
end

%% Moments thresholding algorithm
% requires input image to be of uint8
% returns uint8 threshold
function level = momentsThresh (y)
  n = numel (y) - 1;

  % The threshold is chosen such that partial_sumA(y,t)/partial_sumA(y,n)
  % is closest to x0.
  sumY = sum (y);
  Avec = cumsum (y) / sumY;

  sumB = partial_sumB (y,n);
  sumC = partial_sumC (y,n);
  sumD = partial_sumD (y,n);
  % The following finds x0.
  x2 = (sumB*sumC - sumY*sumD) / (sumY*sumC - sumB^2);
  x1 = (sumB*sumD - sumC^2) / (sumY*sumC - sumB^2);
  x0 = .5 - (sumB/sumY + x2/2) / sqrt (x2^2 - 4*x1);

  % And finally the threshold
  [~, ind] = min (abs (Avec-x0));
  level = ind-1;
end

%% Max Entropy
function T = maxentropy(y)
  n = numel (y) - 1;

  warning ('off', 'Octave:divide-by-zero', 'local');

  % The threshold is chosen such that the following expression is minimized.
  sumY = sum (y);
  negY = negativeE (y, n);
  for j = 0:n
    sumA = partial_sumA (y, j);
    negE = negativeE (y, j);
    sum_diff = sumY - sumA;
    vec(j+1) = negE/sumA - log10 (sumA) + (negY-negE)/(sum_diff) - log10 (sum_diff);
  end

  [~,ind] = min (vec);
  T = ind-1;
end

%% Maximum likelihood
function Tout = maxlikelihood (y)

  n = numel (y) - 1;

  % initial estimate for the threshold is found with the minimum algorithm
  T = minimum (y);

  sumY = sum (y);

  sumB  = partial_sumB (y, n);
  sumC  = partial_sumC (y, n);

  sumAT = partial_sumA (y, T);
  sumBT = partial_sumB (y, T);
  sumCT = partial_sumC (y, T);
  % initial values for the statistics
  mu     = sumBT / sumAT;
  nu     = (sumB - sumBT) / (sumY - sumAT);
  p      = sumAT / sumY;
  q      = (sumY - sumAT) / sumY;
  sigma2 = sumCT / sumAT - mu^2;
  tau2   = (sumC - sumCT) / (sumY - sumAT) - nu^2;

  % Return if sigma2 or tau2 are zero, to avoid division by zero
  if (sigma2 == 0 || tau2 == 0)
    Tout = T;
    return
  end

  not_done = true;
  while not_done
    % we store the previous values for comparison at the end (we will stop when
    % they stop changing)
    mu_prev     = mu;
    nu_prev     = nu;
    p_prev      = p;
    q_prev      = q;
    sigma2_prev = nu;
    tau2_prev   = nu;

    for i = 0:n
      phi(i+1) = p/sqrt((sigma2)) * exp(-((i-mu)^2) / (2*sigma2)) / ...
                (p/sqrt(sigma2)   * exp(-((i-mu)^2) / (2*sigma2)) + ...
                (q/sqrt(tau2))    * exp(-((i-nu)^2) / (2*tau2)));
    end
    ind   = 0:n;
    gamma = 1-phi;
    F     = phi*y';
    G     = gamma*y';

    mu      = ind.*phi*y'/F;
    nu      = ind.*gamma*y'/G;
    p       = F / sumY;
    q       = G / sumY;
    sigma2  = ind.^2.*phi*y'/F - mu^2;
    tau2    = ind.^2.*gamma*y'/G - nu^2;

    not_done =  (abs (mu - mu_prev) <= eps || abs (nu - nu_prev)     <= eps || ...
         abs (p  - p_prev)          <= eps || abs (q - q_prev)       <= eps || ...
         abs (sigma2 - sigma2_prev) <= eps || abs (tau2 - tau2_prev) <= eps);
  end
  % the terms of the quadratic equation to be solved
  w0 = 1/sigma2-1/tau2;
  w1 = mu/sigma2-nu/tau2;
  w2 = mu^2/sigma2 - nu^2/tau2 + log10((sigma2*q^2)/(tau2*p^2));

  % If the threshold would be imaginary, return with threshold set to zero
  sqterm = w1^2-w0*w2;
  if (sqterm < 0)
    Tout{1} = 0;
    return
  end

  % The threshold is the integer part of the solution of the quadratic equation
  Tout = floor((w1+sqrt(sqterm))/w0);
end

%% Minimum thresholding
function T = minimum(y)
  n = numel (y) - 1;

  % Smooth the histogram by iterative three point mean filtering.
  iter = 0;
  while ~bimodtest(y)
    h = ones(1,3)/3;
    y = conv2(y,h,'same');
    iter = iter+1;
    % If the histogram turns out not to be bimodal, set T to zero.
    if iter > 10000;
      T{1} = 0;
      return
    end
  end

  peakfound = false;
  for k = 2:n
    if y(k-1) < y(k) && y(k+1) < y(k)
      peakfound = true;
    end
    if peakfound && y(k-1) >= y(k) && y(k+1) >= y(k)
      T{1} = k-1;
      return
    end
  end
end

%% Minimum error
function [Tout] = minerror_iter (y, T)
  n = numel (y) - 1;

  Tprev = NaN;

  warning ('off', 'Octave:divide-by-zero', 'local');

  sumA = partial_sumA (y, n);
  sumB = partial_sumB (y, n);
  sumC = partial_sumC (y, n);
  while T ~= Tprev
    % Calculate some statistics.
    sumAT = partial_sumA (y, T);
    sumBT = partial_sumB (y, T);
    sumCT = partial_sumC (y, T);
    sumAdiff = sumA - sumAT;

    mu = sumBT/sumAT;
    nu = (sumB-sumBT)/(sumAdiff);
    p = sumAT/sumA;
    q = (sumAdiff) / sumA;
    sigma2 = sumCT/sumAT-mu^2;
    tau2 = (sumC-sumCT) / (sumAdiff) - nu^2;

    % The terms of the quadratic equation to be solved.
    w0 = 1/sigma2-1/tau2;
    w1 = mu/sigma2-nu/tau2;
    w2 = mu^2/sigma2 - nu^2/tau2 + log10((sigma2*q^2)/(tau2*p^2));

    % If the next threshold would be imaginary, return with the current one.
    sqterm = w1^2-w0*w2;
    if sqterm < 0
      warning('MINERROR:NaN','Warning: th_minerror_iter did not converge.')
      break
    end

    % The updated threshold is the integer part of the solution of the
    % quadratic equation.
    Tprev = T;
    T = floor((w1+sqrt(sqterm))/w0);

    % If the threshold turns out to be NaN, return with the previous threshold.
    if isnan(T)
      warning('MINERROR:NaN','Warning: th_minerror_iter did not converge.')
      T = Tprev;
    end
  end
  Tout = T;
end


%% Percentile
% The threshold is chosen such that 50% (in case of p = 0.5) of
% pixels lie in each category.
function [T] = percentile (y)
  p = 0.5;
  Avec = cumsum (y) / sum (y);
  [~, ind] = min (abs (Avec - p));
  T{1} = ind -1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%% Test if a histogram is bimodal.
function b = bimodtest(y)
  len = length(y);
  b = false;
  modes = 0;

  % Count the number of modes of the histogram in a loop. If the number
  % exceeds 2, return with boolean return value false.
  for k = 2:len-1
    if y(k-1) < y(k) && y(k+1) < y(k)
      modes = modes+1;
      if modes > 2
        return
      end
    end
  end

  % The number of modes could be less than two here
  if modes == 2
    b = true;
  end
end

%% Find the local maxima of a vector using a three point neighborhood.
function y = flocmax(x)
%  y    binary vector with maxima of x marked as ones

  len = length(x);
  y = zeros(1,len);

  for k = 2:len-1
    [~,ind] = max(x(k-1:k+1));
    if ind == 2
      y(k) = 1;
    end
  end
end

%% Calculate the balance measure of the histogram around a histogram index.
function E = hbalance(y,ind)
%  y    histogram
%  ind  index about which balance is calculated
%
% Out:
%  E    balance measure
%
% References:
%
% A. Rosenfeld and P. De La Torre, 'Histogram concavity analysis as an aid
% in threshold selection,' IEEE Transactions on Systems, Man, and
% Cybernetics, vol. 13, pp. 231-235, 1983.
%
% P. K. Sahoo, S. Soltani, and A. K. C. Wong, 'A survey of thresholding
% techniques,' Computer Vision, Graphics, and Image Processing, vol. 41,
% pp. 233-260, 1988.

  n = length(y)-1;
  E = partial_sumA(y,ind)*(partial_sumA(y,n)-partial_sumA(y,ind));
end

% Find the convex hull of a histogram.
function H = hconvhull(h)
  % In:
  %  h    histogram
  %
  % Out:
  %  H    convex hull of histogram
  %
  % References:
  %
  % A. Rosenfeld and P. De La Torre, 'Histogram concavity analysis as an aid
  % in threshold selection,' IEEE Transactions on Systems, Man, and
  % Cybernetics, vol. 13, pp. 231-235, 1983.

  len = length(h);
  K(1) = 1;
  k = 1;

  % The vector K gives the locations of the vertices of the convex hull.
  while K(k)~=len

    theta = zeros(1,len-K(k));
    for i = K(k)+1:len
      x = i-K(k);
      y = h(i)-h(K(k));
      theta(i-K(k)) = atan2(y,x);
    end

    maximum = max(theta);
    maxloc = find(theta==maximum);
    k = k+1;
    K(k) = maxloc(end)+K(k-1);

  end

  % Form the convex hull.
  H = zeros(1,len);
  for i = 2:length(K)
    H(K(i-1):K(i)) = h(K(i-1))+(h(K(i))-h(K(i-1)))/(K(i)-K(i-1))*(0:K(i)-K(i-1));
  end
end

% Entropy function. Note that the function returns the negative of entropy.
function x = negativeE(y,j)
  % used by the maxentropy method only
  y = y(1:j+1);
  y = y(y~=0);
  x = sum(y.*log10(y));
end

