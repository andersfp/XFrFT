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
% Number of pixels (must be even)
m = 500;

% Field-of-view
dx = 10e-3;

% Make coordinates
x0 = ((-m/2):(m/2 - 1))/m*dx;
y0 = x0.';
[X0,Y0] = meshgrid(x0,y0);

% Make object
E0 = sqrt(exp(-X0.^2./(2*(400e-6./(2*sqrt(2*log(2)))).^2)).*exp(-Y0.^2./(2*(200e-6./(2*sqrt(2*log(2)))).^2)));


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
d1 = 0;
n = 101;
d2 = linspace(fN-1e-3,fN+1e-3,n);

% Calculate analytical apertures
[sigma_D,sigma_a,sigV,gamma,sigma_p] = CRL_Parameters_2(N,R,mu,f,phi,d1);

% Calculate the effective vignetting width
sigma_v = Vignetting(R,N,mu,d1,T,f,lambda,sigma_p);

% Make an array of the focal lengths
F = f*ones(N,1);

% Set the object plane curvature and scaling parameter
R0 = Inf;
s0 = dx/sqrt(m);

% Calculate the propagation parameters
a = zeros(N+1,n);
Rm = a;
Rp = a;
sm = a;
sp = a;
gm = a;
gp = a;
for i = 1:n
    D = [d1 + T/2;T*ones(N-1,1);T/2 + d2(i)];
    [a(:,i),Rm(:,i),Rp(:,i),sm(:,i),sp(:,i),gm(:,i),gp(:,i)] = FrFT_parameters(D,F,lambda,R0,s0);
end

% Check that the product of the g-parameters is in the range 0<=gm*gp<=1
if any(any(0 > gm.*gp)) || any(any(gm.*gp > 1))
    warning('Product gm*gp is outside limits.');
end


%% Calculate the wave propagation
% Apply the vignetting
E1 = E0.*sqrt(exp(-(X0.^2 + Y0.^2)./(2*sigma_v.^2)));

% Propopagate to the end of the CRL
[E1,X1,Y1] = propFrFT2(E1,X0,Y0,Inf,Inf,sm(1,1),sp(N,1),sum(a(1:N,1)),lambda,sum(D(1:N)),method);

% Apply the effective pupil
E1 = E1.*sqrt(exp(-(X1.^2 + Y1.^2)./(2*sigma_p.^2)));

% Simulate each propagation distance
E2 = zeros(m,m,n);
X2 = E2;
Y2 = E2;
for i = 1:n
    [E2(:,:,i),X2(:,:,i),Y2(:,:,i)] = propFrFT2(E1,X1,Y1,Inf,Inf,sm(end,i),sp(end,i),a(end,i),lambda,T/2 + d2(i),method);
    fprintf('.');
end
fprintf('\n');

% Calculate intensities
I0 = abs(E0).^2;
I1 = abs(E1).^2;
I2 = abs(E2).^2;

% Interpolate the intensity
x = ((-m/2):(m/2 - 1))/m*5e-7;
y = x.';
[X,Y] = meshgrid(x,y);
I = zeros(length(y),length(x),n);
for i = 1:n
    I(:,:,i) = interp2(X2(:,:,i),Y2(:,:,i),I2(:,:,i),X,Y,'spline',0);
    fprintf('.');
end
fprintf('\n');


%% Plot the results
% Plot the incoming beam
figure;
imagesc(1e3*x0,1e3*y0,I0);
axis equal tight;
title('Incoming beam intensity');
xlabel('x [mm]');
ylabel('y [mm]');

% Slice plots of the focal point
figure;
slice(1e9*x,1e9*y,1e6*(d2 - fN),I,0,0,[0,200]);
shading flat;
axis equal tight;
xlabel('x [nm]');
ylabel('y [nm]');
zlabel('Propagation distance from focal point [\mum]');

% Isosurface plot of the focal point
cm = winter(11);
figure;
p1 = patch(isosurface(1e9*x,1e9*y,1e6*(d2 - fN),I,0.2*max(I(:))));
p1.FaceColor = cm(3,:);
p1.EdgeColor = 'none';
p1.FaceAlpha = 0.2;
hold on;
p2 = patch(isosurface(1e9*x,1e9*y,1e6*(d2 - fN),I,0.5*max(I(:))));
p2.FaceColor = cm(6,:);
p2.EdgeColor = 'none';
p2.FaceAlpha = 0.5;
p3 = patch(isosurface(1e9*x,1e9*y,1e6*(d2 - fN),I,0.8*max(I(:))));
p3.FaceColor = cm(9,:);
p3.EdgeColor = 'none';
p3.FaceAlpha = 0.8;
p4 = patch(isosurface(1e9*x,1e9*y,1e6*(d2 - fN),I,0.1*max(I(:))));
p4.FaceColor = cm(2,:);
p4.EdgeColor = 'none';
p4.FaceAlpha = 0.1;
axis equal;
camlight;
lighting gouraud;
xlabel('x [nm]');
ylabel('y [nm]');
zlabel('Propagation distance from focal point [\mum]');


