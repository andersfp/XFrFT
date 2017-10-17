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
% Number of pixels on each side (must be even)
m = 1000;

% Field-of-view
dx = 10e-6;

% Make coordinates
x0 = ((-m/2):(m/2 - 1))/m*dx;
y0 = x0.';
[X0,Y0] = meshgrid(x0,y0);

% Make object
E0 = zeros(m,m);
E0(m/2+1,m/2+1) = 1; % Center pixel at (0,0) is set to 1

% Define the slit sizes to test
n = 40;
sl = logspace(-4.5,-2.65,n).';


%% Set the simulation details
% Get the wavelength from the X-ray energy
E = 17e3;
lambda = 1e-10*12398.42./E;

% Get material properties (delta and mu) at this energy
[delta,mu] = Be_Prop(E);

% Set CRL parameters
R = 50e-6; % Radius of curvature at apex
T = 1.6e-3; % Distance between each lens element
Tweb = 2e-6; % Minimum Be thickness at the optical axis
N = 69; % Number of lenses

% Calculate CRL focal lengths
[f,phi,fN] = CRL_Parameters_1(R,T,N,delta);

% Set object and image plane positions (see CRL paper)
M = 10;
d1 = fN.*(1 + 1./(M.*cos(N.*phi)));
d2 = fN.*(1 + M./cos(N.*phi));

% Calculate analytical apertures
[sigma_D,sigma_a,sigV,gamma,sigma_p] = CRL_Parameters_2(N,R,mu,f,phi,d1);

% Calculate the effective vignetting width
sigma_v = Vignetting(R,N,mu,d1,T,f,lambda,sigma_p);

% Make an array of the focal lengths
F = [f*ones(N,1);Inf]; % We need to access the BFP, and we do that by adding a fictional lens with an infinite focal length

% Make an array of distances in between each lens (plus object-lens and lens-detector)
D = [d1 + T/2;T*ones(N-1,1);T/2 + fN;d2 - fN]; % We need to access the BFP, and we do that by adding a fictional lens at the location of the BFP

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
% Define function for attenuation (CRL paper, sqrt(...) due to amplitude)
att = @(x,y) sqrt(exp(-mu*(x.^2 + y.^2)./R));

% Propagate to the BFP, where the slit is located
E1 = E0.*sqrt(exp(-(X0.^2 + Y0.^2)./(2*sigma_v.^2))); % Apply the vignetting (sqrt(...) due to amplitude)
[E1,X1,Y1] = propFrFT2(E1,X0,Y0,Inf,Inf,sm(1),sp(end-2),sum(a(1:end-2)),lambda,sum(D(1:end-2)),method); % Propagate to the end of the CRL in 1 transform by summing the a-parameter
E1 = E1.*sqrt(exp(-(X1.^2 + Y1.^2)./(2*sigma_p.^2))); % Apply the effective pupil (sqrt(...) due to amplitude)
[E2,X2,Y2] = propFrFT2(E1,X1,Y1,Inf,Inf,sm(end-1),sp(end-1),a(end-1),lambda,D(end-1),method); % Propagate to the BFP

% Apply the slit and propagate to the BFP
ED = zeros(m,m,n);
tic;
for i = 1:n
    Ed = E2.*rect(X2/sl(i)).*rect(Y2/sl(i)); % Apply the slit
    [ED(:,:,i),XD,YD] = propFrFT2(Ed,X2,Y2,Inf,Inf,sm(end),sp(end),a(end),lambda,D(end),method); % Propagation to detector
    fprintf('.');
end
fprintf('\n');
toc;


%% Analyze the PSF
% Get the axes
x1 = X1(1,:);
y1 = Y1(:,1);
x2 = X2(1,:);
y2 = Y2(:,1);
xD = XD(1,:);
yD = YD(:,1);

% Calibrate the slit sizes (due to pixel quantization)
sl = sum(rect(repmat(x2,n,1)./repmat(sl,1,m)),2).*mean(diff(x2));

% Calculate intensities
ID = abs(ED).^2;

% Get vertical line profiles through the optical axis
L = squeeze(ID(:,m/2+1,:));

% Normalize the line profiles
Ln = L./repmat(max(L),m,1);

% Geometrical optics PSF widths
sigma_crl = repmat(lambda./(4*pi*sigma_a),n,1);
sigma_slit = 0.3645*lambda*fN./(sl*cos(N*phi));
sigma_tot = sqrt(sigma_crl.^2 + sigma_slit.^2);

% Fit the simulated PSFs with a Gaussian
fun = @(sigma,x) exp(-x.^2./(2*sigma.^2));
sigma_frft = zeros(n,1);
for i = 1:n
    ff = fit(1e8*yD/M,Ln(:,i),fun,'StartPoint',1e8*sigma_tot(i));
    sigma_frft(i) = ff.sigma/1e8;
end

% Plot the results
figure;
loglog(1e3*sl,1e9*[sigma_crl,sigma_slit,sigma_tot,sigma_frft],'LineWidth',2);
set(gca,'XLim',[3e-2 3],'YLim',[4 300]);
title('PSF width');
xlabel('Slit size [mm]');
ylabel('PSF rms width [nm]');
legend('\sigma_{CRL}','\sigma_{slit}','\sigma_{tot}','\sigma_{FrFT}');



