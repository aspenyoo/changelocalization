function gabor = makeGabor(imSize, contrast, theta, wavelength, phaseShift)
% 
% outputs matrix with values 0 to 255
% 
% ===== INPUT VARIABLES ======
% IMSIZE: scalar. image size in pixels imSize x imSize
% CONTRAST: Michelson luminance contrast. value 0 to 1. default 1
% ORIENTATION: 0 to 180 deg. default 0
% WAVELENGTH: how many pixels per cycle. default IMSIZE/5
% PHASESHIFT: proportion phase shift rightward. 0 to 1. default 0

if nargin < 2; contrast = 1; end
if nargin < 3; theta = 0; end
if nargin < 4; wavelength = imSize/5; end % default freq of 5?
if nargin < 5; phaseShift = 0; end

% clear all
% contrast = 1;                           % contrast: ranges from 0 to 1
% imSize = 100;                           % image size: n X n
% wavelength = 30;                             % wavelength (number of pixels per cycle)
% theta = 0;                             % grating orientation in degrees
% phaseShift = 0;                            % phase (0 -> 1)

sigma = 0.15;                             % guassian SD as proportion of imSize
trim = .005;                            % trim off gaussian values smaller than this
spatialFreq = imSize/wavelength;                    % compute frequency from wavelength
thetaRad = deg2rad(theta);        % convert theta (orientation) to radians
phaseRad = 2*pi*phaseShift;             % convert to radians: 0 -> 2*pi

X = linspace(-.5,.5,imSize);
[Xm, Ym] = meshgrid(X, X);              
Xt = Xm * sin(thetaRad);                % compute proportion of Xm for given orientation
Yt = Ym * cos(thetaRad);                % compute proportion of Ym for given orientation
XYt = Xt + Yt;                          % sum X and Y components
XYf = XYt * spatialFreq * 2*pi;                % convert to radians and scale by frequency
grating = sin( XYf + phaseRad);                   % make 2D sinewave

gauss = exp(-((Xm.^2)+(Ym.^2))./(2* sigma^2));      % make 2D gaussian
gauss(gauss < trim) = 0;                 % trim around edges

gabor = round(255.*(contrast.*(grating .* gauss)+1)./2);       % element wise multiplication.

% imagesc( gabor, [0 1] );    colormap gray(256);    