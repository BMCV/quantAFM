%----------------------------------------------
% This is a config file for the ANuFE project
% All user specific settings should be set here
% Care that only images with 'quadratic' pixels are supported, i.e. they
% have the same length in x- and y- direction
% All pictures processed in one run have to have the same x-resolution and
% scansize
%----------------------------------------------

% Main picture directory, all pictures should be in the respective
% subfolders with regard to the current operating system
workingDirWindows = '..\denoised_imgs';
workingDirLinux = '../denoised_imgs';
% Respective subfolder whose pictures are to be analyzed with regard to the
% current operating system
currentImageDirWindows = '..\denoised_imgs\p_Wildtyp\*.tif';
currentImageDirLinux = '../denoised_imgs/p_Wildtyp/';

% Set to 1 if the main DNA backbone should be reconstructed after thinning 
% (i.e. recovery of 'lost' pixels ,0 else
recoverBackbone = 0; %Default: 1

% set min and max length for recognized DNA length without histones in nm
minLength_free = 170; % Default: 180
maxLength_free = 230; % Default: 250

% set min and max length for recognized DNA length with attached histones in nm
minLength_bound = 90; % Default: 100
maxLength_bound = 210; % Default: 200

% set min and max bounds for radius of nuclei in nm
minRadius = 9.375; % Default: 9.375
maxRadius = 14.0625; % Default: 14.0625

% set length of the picture scanline in x-Direction, i.e. from left to
% right in �m
scansize = 3;

% give the x-resolution of the images
xResolution = 1280;

% only output objects for which all calculations ended successfully if 1, 0
% else
onlyValid = 1; % Default: 1

% 1 if you want to recalculate the global threshold presets, 0 else
% Set this only to 0 if you already processed your dataset
setMeanThreshold = 1; % Default: 1

% Use multiple CPU cores, recommended for big datasets
% 1 for multiple cores, 0 for single core
parallel = 0; % Default: 0

% 1 when using GPU, 0 else
% only use if you have a CUDA compatible NVIDIA graphics card
% use with caution and disable if it does not work, there was no test which
% CUDA versions are supported
gpu = 0; % Default: 0

%----------------------------------------------
% WARNING: All of the following is only for experienced users/developers

% The threshold algorithms used in the three step thresholding algorithm,
% more information in the project report, commonly used are 'otsu' and
% 'moments', other Matlab algorithms did not work for our dataset
threshAlgorithm1 = 'otsu'; % Default: 'otsu'
threshAlgorithm2 = 'otsu'; % Default: 'otsu'
% Upper bound for background noise, every value that is <=
% backgroundThreshold is considered background
backgroundThreshold = 95; % Default: 95


% Additionally, you can specify the used algorithm for length determination
% lengthAlgo = 'C'; % default arbitrary
lengthAlgo = 'D';


% determine whether the output is "real-valued" (in nm), or "pixel-valued"
realValue = 1; % default is 1 = real value
% realValue = 0; % this is pixel values
%----------------------------------------------