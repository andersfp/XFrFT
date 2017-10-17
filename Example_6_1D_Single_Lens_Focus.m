% Initialization (optional)
clear;
close all;
clc;


%% Define the object
% Number of pixels (must be even)
m = 10000;

% Field-of-view
dx = 100e-6;

% Make coordinates
x0 = ((-m/2):(m/2 - 1)).'/m*dx;

% Make object
w = 10e-6; % Lens pupil width
E0 = rect(x0/w);


%% Set the simulation details
% Get the wavelength from the X-ray energy
E = 17e3;
lambda = 1e-10*12398.42./E;

% Determine propagation distance steps after lens
n = 1000;
z = linspace(0,0.03,n);

% Set the focal length
F = 0.02;

% Set the object plane curvature and scaling parameter
R0 = Inf;
s0 = dx/sqrt(m);

% Calculate the propagation parameters
a = zeros(2,n);
Rm = a;
Rp = a;
sm = a;
sp = a;
gm = a;
gp = a;
for i = 1:n
    [a(:,i),Rm(:,i),Rp(:,i),sm(:,i),sp(:,i),gm(:,i),gp(:,i)] = FrFT_parameters([0;z(i)],F,lambda,R0,s0);
end

% Check that the product of the g-parameters is in the range 0<=gm*gp<=1
if any(any(0 > gm.*gp)) || any(any(gm.*gp > 1))
    warning('Product gm*gp is outside limits.');
end


%% Calculate the wave propagation
% Simulate each propagation distance
E1 = zeros(m,n);
x1 = E1;
for i = 1:n
    [E1(:,i),x1(:,i)] = propFrFT1(E0,x0,Rm(1,i),Rp(2,i),sm(1,i),sp(2,i),sum(a(:,i)),lambda,z(i));
    if mod(i,10) == 0
        fprintf('.');
    end
end
fprintf('\n');

% Calculate intensities
I0 = abs(E0).^2;
I1 = abs(E1).^2;

% Interpolate the intensity
xmin = min(min(x1));
xmax = max(max(x1));
dx1 = min(min(abs(diff(x1))));
x = (xmin:dx1:xmax).';
I = zeros(length(x),n);
for i = 1:n
    I(:,i) = interp1(x1(:,i),I1(:,i),x,'spline',0);
    if mod(i,10) == 0
        fprintf('.');
    end
end
fprintf('\n');


%% Plot the results
% Plot the gradual propagation
figure;
imagesc(1e2*z,1e6*x,I);
xlabel('Propagation distance [cm]');
ylabel('x [\mum]');

% Plot the gradual propagation (zoom)
figure;
imagesc(1e2*z,1e6*x,I);
set(gca,'YLim',[-10 10]);
xlabel('Propagation distance [cm]');
ylabel('x [\mum]');

% Plot the transform order parameter
figure;
plot(1e2*z,sum(a));
xlabel('Propagation distance [cm]');
ylabel('Transform order');



