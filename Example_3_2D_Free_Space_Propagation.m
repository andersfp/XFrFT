% Initialization (optional)
clear;
close all;
clc;


%% Decide how to calculate the FrFT
% Define method
% method = 'gpu'; % Fastest implementation, but requires an Nvidia GPU.
method = 'vec'; % Fastest implementation on CPU (for not too large images), based on vectorized array operations.
% method = 'par'; % Parallel for-loop implementation on CPU, can be faster than 'vec' for large images, and uses less memory than 'vec' requires parallel toolbox.
% method = 'for'; % Non-parallel for-loop on CPU, the slowest implementation, but has the lowest memory consumption and does not require parallel toolbox.


%% Define the object
% Get a test image
I0 = double(imread('Matterhorn.png')); % Picture from pixabay.com

% Zeropad the image (Best to have 0 intensity close to the edges)
I0 = padarray(I0,[50 50],0,'both');

% Get the amplitude from the intensity
E0 = sqrt(I0);

% Number of pixels (must be even)
m = size(E0,1);

% Field-of-view
dx = 100e-6;

% Make coordinates
x0 = ((-m/2):(m/2 - 1))/m*dx;
y0 = x0.';
[X0,Y0] = meshgrid(x0,y0);


%% Set the simulation details
% Get the wavelength from the X-ray energy
E = 17e3;
lambda = 1e-10*12398.42./E;

% Focal length is an empty array, as we do not have any lenses
F = [];

% Propagation distance
D = 0.2;

% Set the object plane curvature and scaling parameter
R0 = Inf;
s0 = dx/sqrt(m);

% Calculate the propagation parameters
[a,Rm,Rp,sm,sp,gm,gp] = FrFT_parameters(D,F,lambda,R0,s0);

% Check that the product of the g-parameters is in the range 0<=gm*gp<=1
if any(0 > gm.*gp) || any(gm.*gp > 1)
    warning('Product gm*gp is outside limits.');
end


%% Calculate the wave propagation
% Free space propagation
[E1,X1,Y1] = propFrFT2(E0,X0,Y0,Rm,Rp,sm,sp,a,lambda,D,method);

% Calculate intensities
I1 = abs(E1).^2;


%% Plot the results
% Get the axes
x1 = X1(1,:);
y1 = Y1(:,1);

% Plot the initial image
figure;
imagesc(1e6*x0,1e6*y0,I0);
axis equal tight;
colormap gray;
title('Object');
xlabel('x [\mum]');
ylabel('y [\mum]');

% Plot the detector image
figure;
imagesc(1e6*x1,1e6*y1,I1);
axis equal tight;
colormap gray;
title('Detector');
xlabel('x [\mum]');
ylabel('y [\mum]');


