% Initialization (optional)
clear;
close all;
clc;


%% Define the object
% Number of pixels on each side (must be even)
m = 1000;

% Field-of-view
dx = 10e-6;

% Make coordinates
x0 = ((-m/2):(m/2 - 1)).'/m*dx;

% Make object
E0 = zeros(m,1);
E0(m/2+1) = 1;


%% Set the simulation details
% Set the central energy
E = 17e3;

% Set the highest energy RMS bandwidth
dEmax = 1e-2*E;

% Set the energy steps
dE = 1;

% Generate energy and wavelength samples
Es = ((-3*dEmax):1:(3*dEmax)) + E;
lambda = 1e-10*12398.42./Es;
n = length(Es);
n0 = round(n/2);

% Define the energy RMS widths
w = [1e-10 (1:100)*1e-4];
w = permute(w,[3 1 2]);
sigma_E = w.*E;
k = length(w);

% Get material properties (delta and mu) at this energy
[delta,mu] = Be_Prop(Es);

% Set CRL parameters
R = 50e-6; % Radius of curvature at apex
T = 1.6e-3; % Distance between each lens element
Tweb = 2e-6; % Minimum Be thickness at the optical axis
N = 69; % Number of lenses

% Calculate CRL focal lengths
[f,phi,fN] = CRL_Parameters_1(R,T,N,delta);

% Set object and image plane positions (see CRL paper)
M = 10;
d1 = fN(n0).*(1 + 1./(M.*cos(N.*phi(n0))));
d2 = fN(n0).*(1 + M./cos(N.*phi(n0)));

% Calculate analytical apertures
[sigma_D,sigma_a,sigV,gamma,sigma_p] = CRL_Parameters_2(N,R,mu,f,phi,d1);

% Calculate the effective vignetting width
sigma_v = zeros(1,n);
for i = 1:n
    sigma_v(i) = Vignetting(R,N,mu(i),d1,T,f(i),lambda(i),sigma_p(i));
    if mod(i,10) == 0
        fprintf('.');
    end
end
fprintf('\n');

% Make an array of distances in between each lens (plus object-lens and lens-detector)
D = [d1 + T/2;T*ones(N-1,1);T/2 + d2];

% Set the object plane curvature and scaling parameter
R0 = Inf;
s0 = dx/sqrt(m);


%% Calculate the wave propagation
% Perform the propagation
E2 = zeros(m,n);
x2 = E2;
for i = 1:n
    F = f(i)*ones(N,1);
    [a,Rm,Rp,sm,sp,gm,gp] = FrFT_parameters(D,F,lambda(i),R0,s0);
    if any(0 > gm.*gp) || any(gm.*gp > 1)
        warning('Product gm*gp is outside limits.');
    end
    E1 = E0.*sqrt(exp(-x0.^2./(2*sigma_v(i).^2)));
    [E1,x1] = propFrFT1(E1,x0,Inf,Inf,sm(1),sp(end-1),sum(a(1:end-1)),lambda(i),sum(D(1:end-1)));
    E1 = E1.*sqrt(exp(-x1.^2./(2*sigma_p(i).^2)));
    [E2(:,i),x2(:,i)] = propFrFT1(E1,x1,Inf,Inf,sm(end),sp(end),a(end),lambda(i),D(end));
    if mod(i,10) == 0
        fprintf('.');
    end
end
fprintf('\n');


%% Process the raw results
% Get the intensities
I0 = abs(E0).^2;
I2 = abs(E2).^2;

% Plot the detector plane
figure;
pcolor(1e-3*Es,1e6*x2,I2);
shading flat;
title('Detector intensity');
xlabel('Energy [keV]');
ylabel('x [\mum]');

% Interpolate the data onto a common grid
x = x2(:,n0)/M;
I = zeros(m,n);
for i = 1:n
    I(:,i) = interp1(x2(:,i),I2(:,i),x*M,'spline',0);
    if mod(i,10) == 0
        fprintf('.');
    end
end
fprintf('\n');

% Plot the new detector plane
figure;
imagesc(1e-3*Es,1e6*x,I);
title('Detector intensity');
xlabel('Energy [keV]');
ylabel('x [\mum]');

% Generate intensity weighting factor based on the energy bandwidths
S = exp(-(repmat(Es,1,1,k) - E).^2./(2.*repmat(sigma_E,1,n,1).^2));

% Generate total intensities
psft = squeeze(sum(repmat(I,1,1,k).*repmat(S,m,1,1),2));

% Plot the total PSFs
figure;
imagesc(1e2*squeeze(w),1e6*x,psft);
title('Simulated PSF vs. bandwidth');
xlabel('RMS bandwidth [%]');
ylabel('x [\mum]');

% Calculate normalized PSFs
%psfn = psft./repmat(sum(psft),m,1); % Area normalized
psfn = psft./repmat(max(psft),m,1); % Maximum normalized

% Plot the normalized intensities
figure;
imagesc(1e2*squeeze(w),1e6*x,psft);
title('Simulated normalized PSF vs. bandwidth');
xlabel('RMS bandwidth [%]');
ylabel('x [\mum]');

% Calculate the RMS width of the PSFs
sigma_psf = sqrt(sum(psft.*repmat(x,1,k).^2)./sum(psft));

% Plot the RMS width
figure;
plot(1e2*squeeze(w),1e9*sigma_psf);
title('RMS PSF width');
xlabel('RMS bandwidth [%]');
ylabel('RMS width [nm]');

% Plot selected PSFs
figure;
plot(1e9*x,psfn(:,[1 2 11 101]),'LineWidth',2);
set(gca,'XLim',[-500 500]);
title('Selected PSFs');
xlabel('Position [nm]');
ylabel('Normalized intensity');
legend('\DeltaE/E = 0','\DeltaE/E = 10^{-4}','\DeltaE/E = 10^{-3}','\DeltaE/E = 10^{-2}');


%% Comparison to geometrical optics
% Calculate parameters
dch = N*phi(n0)*(d1*d2/(f(n0)*phi(n0)) - f(n0)*phi(n0))*cos(N*phi(n0)) + (d1*d2/(f(n0)*phi(n0)) + f(n0)*phi(n0) + N*phi(n0)*(d1 + d2))*sin(N*phi(n0));

% Calculate the chromatic PSF
psfgt = squeeze(exp(-abs(repmat(x*M,1,1,k))./(sigma_a(n0).*repmat(w,m,1,1).*dch)));

% Normalize the PSF
%psfgn = psfgt./repmat(sum(psfgt),m,1); % Area normalized
psfgn = psfgt./repmat(max(psfgt),m,1); % Maximum normalized

% Plot the PSFs versus bandwidth
figure;
subplot(1,2,1);
imagesc(100*squeeze(w),1e9*x,psft/max(max(psft)),[0 1]);
hold on;
[C,h] = contour(100*squeeze(w),1e9*x,psfn,[0.25 0.5 0.75],'Color','black');
clabel(C,h,'FontSize',12,'LabelSpacing',170);
xlabel('Relative bandwidth [%]');
ylabel('Position [nm]');
set(gca,'YLim',[-500 500]);
subplot(1,2,2);
imagesc(100*squeeze(w),1e9*x,psfgt/max(max(psfgt)),[0 1]);
hold on;
[C,h] = contour(100*squeeze(w),1e9*x,psfgn,[0.25 0.5 0.75],'Color','black');
clabel(C,h,'FontSize',12,'LabelSpacing',110);
xlabel('Relative bandwidth [%]');
ylabel('Position [nm]');
set(gca,'YLim',[-500 500]);



