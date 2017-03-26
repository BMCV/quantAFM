%----------------------------------------------
% This is a config file for the ANuFE project
% All user specific settings should be set here
% Care that only images with 'quadratic' pixels are supported, i.e. they
% have the same length in x- and y- direction
% All pictures processed in one run have to have the same x-resolution and
% scansize. Otherwise alter the file for each batch of images that have the
% same scansize and resolution.
%----------------------------------------------

% Main picture directory, all pictures should be in the respective
% subfolders with regard to the current operating system
workingDirWindows = '..\denoised_imgs';
workingDirLinux = '../denoised_imgs';
% Respective subfolder whose pictures are to be analyzed with regard to the
% current operating system
currentImageDirWindows = '..\denoised_imgs\test\*.tif';
currentImageDirLinux = '../denoised_imgs/p_Wildtyp/';

% Set to 1 if the main DNA backbone should be reconstructed after thinning 
% (i.e. recovery of 'lost' pixels ,0 else
recoverBackbone = 1; %Default: 1

% set min and max length for recognized DNA length without histones in nm
minLength_free = 100; % Default: 180
maxLength_free = 530; % Default: 250

% set min and max length for recognized DNA length with attached histones in nm
minLength_bound = 90; % Default: 100
maxLength_bound = 510; % Default: 200

% set min and max bounds for radius of nuclei in nm
minRadius = 9.375; % Default: 9.375
maxRadius = 14.0625; % Default: 14.0625

% set length of the picture scanline in x-Direction, i.e. from left to
% right in mikrometers
scansize = 10;

% give the x-resolution of the images
xResolution = 4096;

% If the images are known to contain dirt of a certain size, they can be
% excluded via the following min/max values.
% This is also relevant if the size of the DNA is already given, since this
% actively affects the preprocessing steps, while the DNA length does not
% (it will only come into play once the length of the objects is
% determined).
% If the values are set to the default of -1, the size will be determined
% by the factors 0.05/0.85 * xResolution
minSize = -1; % -1 default, else positive integer
% minSize = 100;
maxSize = -1; % -1 default, else positive integer
% maxSize = 900;
% As some advice: The values for 660 bp DNA strands have been 100/900 (see
% report)

% only output objects for which all calculations ended successfully if 1, 0
% else
onlyValid = 0; % Default: 1

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

% output folder
% the output csv files will be placed in the appropriate subfolder. This
% serves the compact view of files, and avoids file cluttering
% default value is "output"
outputDirWin = 'output\';
outputDirLinux = 'output/';

%----------------------------------------------
% WARNING: All of the following settings are only for the experienced user

% The threshold algorithms used in the three step thresholding algorithm,
% more information in the project report, commonly used are 'otsu' and
% 'moments', other Matlab algorithms did not work for our dataset
threshAlgorithm1 = 'otsu'; % Default: 'otsu'
threshAlgorithm2 = 'otsu'; % Default: 'otsu'
% Upper bound for background noise, every value that is <=
% backgroundThreshold is considered background
backgroundThreshold = 25; % Default: 95


% Additionally, you can specify the used algorithm for length determination
lengthAlgo = 'C'; 
% lengthAlgo = 'D';
% Method C works with less variance on the test data, and is generally a
% bit more robust,
% while method D has a more restrictive policy for self-intersections and
% objects that would possibly be valid with method C (advisable for dirty
% and noisy datasets)


% determine whether the output is "real-valued" (in nm), or "pixel-valued"
realValue = 1; % default is 1 = real value
% realValue = 0; % this is pixel values


% export the pixel positions of each DNA object to calculate bending
% stiffness from it. This feature was specifically implemented for the
% Langowski group, DKFZ (German Cancer Research Center).
% The output is of the following format:
% x-position y-position
% for each pixel in the connected (thinned) component of each individual
% DNA strand
exportPixels = 1; % either 0 or 1
exportOnlyValid = 0; % either 0 or 1
exportReal = 1; % either 0 r 1
% exportPixels defines whether or not the bending stiffness is of interest
% at all, and whether the respective files should be output or not.
% exportOnlyValid allows the user to control the output for invalid DNA
% objects. For some cases, even the output of these could be interesting
% exportReal defines the measure of the values. If the provided value is 1,
% the measures will be converted to nanometers/mikrometeres, otherwise
% pixel values.

% the default subfolder is "export".
exportDirWin = 'export\';
exportDirLinux = 'export/';

% Attention! The following value is only experimental!
dkfz = 1;

%----------------------------------------------