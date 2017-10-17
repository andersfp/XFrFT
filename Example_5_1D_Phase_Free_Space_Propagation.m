% Initialization (optional)
clear;
close all;
clc;


%% Define the object
% Number of pixels (must be even)
m = 1000;

% Field-of-view
dx = 100e-6;

% Make coordinates
x0 = ((-m/2):(m/2 - 1)).'/m*dx;

% Make object
w = 50e-6;
E0 = rect(x0/w);


%% Set the simulation details
% Get the wavelength from the X-ray energy
E = 17e3;
lambda = 1e-10*12398.42./E;

% Focal length is an empty array, as we do not have any lenses
F = [];

% Propagation distance
D = 1;

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
[E1,x1] = propFrFT1(E0,x0,Inf,Inf,sm,sp,a,lambda,D);
% We are interested in the phase, and since we do not have a very fine
% sampling grid we do not apply the quadratic phase in the detector plane.
% But we know what the quadratic phase curvature is, so we can add it
% later.

% Calculate intensities (these do not change whether we apply the quadratic
% phase in the detector plane or not)
I0 = abs(E0).^2;
I1 = abs(E1).^2;

% Get the phase in the detector plane (without the quadratic phase)
p1 = unwrap(angle(E1));

% Add a constant phase to set it to 0 at x=0, and add the quadratic phase
% (equation (4b))
p1 = p1 - p1(m/2 + 1) + pi*x1.^2./(lambda*Rp);


%% Plot the results
% Plot the initial and final intensity
figure;
plot(1e6*x0,I0,1e6*x1,I1);
title('Object and detector intensity');
xlabel('x [\mum]');
ylabel('Intensity [a.u.]');

% Plot the phase
figure;
plot(1e6*x1,p1);
title('Detector plane phase');
xlabel('x [\mum]');
ylabel('Phase [rad]');

% Plot the phase (zoom)
figure;
plot(1e6*x1,p1);
set(gca,'XLim',[-50 50],'Ylim',[-pi 3*pi]);
title('Detector plane phase (zoom)');
xlabel('x [\mum]');
ylabel('Phase [rad]');



