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
dx = 1e-6;

% Make coordinates
x0 = ((-m/2):(m/2 - 1))/m*dx;
y0 = x0.';
[X0,Y0] = meshgrid(x0,y0);

% Make object
E0 = zeros(m,m);
E0(m/2+1,m/2+1) = 1; % Center pixel at (0,0) is set to 1


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
F = f*ones(N,1);

% Make an array of distances in between each lens (plus object-lens and lens-detector)
D = [d1 + T/2;T*ones(N-1,1);T/2 + d2];

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

% Propagate in N+1 steps with attenuation
E70 = E0;
X70 = X0;
Y70 = Y0;
tic;
for i = 1:N+1
    [E70,X70,Y70] = propFrFT2(E70,X70,Y70,Inf,Inf,sm(i),sp(i),a(i),lambda,D(i),method); % Important to use Inf for both Rm and Rp, since the transform order already takes the lenses into account
    if i <= N
        E70 = E70.*att(X70,Y70); % Apply the attenuation at each lens (after each of the first N steps, don't apply in the detector plane)
    end
    fprintf('.');
end
fprintf('\n');
toc;

% Propagate in 2 steps
tic;
E1 = E0.*sqrt(exp(-(X0.^2 + Y0.^2)./(2*sigma_v.^2))); % Apply the vignetting (sqrt(...) due to amplitude)
[E1,X1,Y1] = propFrFT2(E1,X0,Y0,Inf,Inf,sm(1),sp(end-1),sum(a(1:end-1)),lambda,sum(D(1:end-1)),method); % Propagate to the end of the CRL in 1 transform by summing the a-parameter
E1 = E1.*sqrt(exp(-(X1.^2 + Y1.^2)./(2*sigma_p.^2))); % Apply the effective pupil (sqrt(...) due to amplitude)
[E2,X2,Y2] = propFrFT2(E1,X1,Y1,Inf,Inf,sm(end),sp(end),a(end),lambda,D(end),method); % Propagate to the detector plane
toc;


%% Plot the results
% Get the axes
x70 = X70(1,:);
y70 = Y70(:,1);
x2 = X2(1,:);
y2 = Y2(:,1);

% Get the intensities
I0 = abs(E0).^2;
I70 = abs(E70).^2;
I2 = abs(E2).^2;

% Get lineplot through the center
l0 = abs(E0(:,m/2+1)).^2;
l70 = abs(E70(:,m/2+1)).^2;
l2 = abs(E2(:,m/2+1)).^2;

% Get phase lineplot through the center
p0 = unwrap(angle(E0(:,m/2+1)));
p70 = unwrap(angle(E70(:,m/2+1)));
p70 = p70 - p70(m/2+1) + pi*y70.^2./(lambda.*Rp(end));
p2 = unwrap(angle(E2(:,m/2+1)));
p2 = p2 - p2(m/2+1) + pi*y2.^2./(lambda.*Rp(end));

% Plot the object plane
figure;
imagesc(x0,y0,I0);
axis equal tight;
title('Object intensity');
xlabel('x [m]');
ylabel('y [m]');

% Plot the detector plane with lens-by-lens simulation
figure;
imagesc(x70,y70,I70);
axis equal tight;
title('Detector intensity (lens-by-lens)');
xlabel('x [m]');
ylabel('y [m]');

% Plot the detector plane with effective vignetting and pupil
figure;
imagesc(x2,y2,I2);
axis equal tight;
title('Detector intensity (effective vignetting and pupil)');
xlabel('x [m]');
ylabel('y [m]');

% Plot the detector plane intensity and phase profiles
figure;
subplot(2,2,1);
plot(x2,l2,'-',x70,l70,'--');
set(gca,'XLim',[-3e-6 3e-6],'YLim',[0 3.5e-10]);
title('Detector intensity');
xlabel('x [m]');
ylabel('Intensity [a.u.]');
legend('Effective aperture','Lens-by-lens');
subplot(2,2,2);
yyaxis left;
plot(x2,l2 - l70);
set(gca,'XLim',[-3e-6 3e-6],'YLim',[-2e-14 2e-14]);
title('Detector intensity error');
xlabel('x [m]');
ylabel('Absolute error [a.u.]');
yyaxis right;
plot(x2,(l2 - l70)./l70);
set(gca,'XLim',[-3e-6 3e-6],'YLim',[-0.01 0.01]);
ylabel('Relative error');
subplot(2,2,3);
plot(x2,p2,'-',x70,p70,'--');
set(gca,'XLim',[-3e-6 3e-6],'YLim',[0 4]);
title('Detector phase');
xlabel('x [m]');
ylabel('Phase [rad]');
legend('Effective aperture','Lens-by-lens');
subplot(2,2,4);
yyaxis left;
plot(x2,p2 - p70);
set(gca,'XLim',[-3e-6 3e-6],'YLim',[-0.5 0.5]);
title('Detector phase error');
xlabel('x [m]');
ylabel('Absolute error [rad]');
yyaxis right;
plot(x2,(p2 - p70)./p70);
set(gca,'XLim',[-3e-6 3e-6],'YLim',[-0.01 0.01]);
ylabel('Relative error');


