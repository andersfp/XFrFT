% Initialization (optional)
clear;
close all;
clc;


%% Define the object
% Number of pixels (must be even)
m = 1000;

% Field-of-view
dx = 10e-3;

% Make coordinates
x0 = ((-m/2):(m/2 - 1)).'/m*dx;

% Make object
E0 = ones(m,1);


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
n = 1001;
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
E1 = E0.*sqrt(exp(-x0.^2./(2*sigma_v.^2)));

% Propopagate to the end of the CRL
[E1,x1] = propFrFT1(E1,x0,Inf,Inf,sm(1,1),sp(N,1),sum(a(1:N,1)),lambda,sum(D(1:N)));

% Apply the effective pupil
E1 = E1.*sqrt(exp(-x1.^2./(2*sigma_p.^2)));

% Simulate each propagation distance
E2 = zeros(m,n);
x2 = E2;
for i = 1:n
    [E2(:,i),x2(:,i)] = propFrFT1(E1,x1,Inf,Inf,sm(end,i),sp(end,i),a(end,i),lambda,0); % Using a propagation distance of 0 to avoid the phase change associated with free space propagation, as this evolves very rapidly.
    if mod(i,10) == 0
        fprintf('.');
    end
end
fprintf('\n');

% Calculate intensities
I0 = abs(E0).^2;
I1 = abs(E1).^2;
I2 = abs(E2).^2;

% Interpolate the intensity
x = ((-10*m/2):(10*m/2 - 1)).'/m*2e-6;
I = zeros(length(x),n);
for i = 1:n
    I(:,i) = interp1(x2(:,i),I2(:,i),x,'spline',0);
    if mod(i,10) == 0
        fprintf('.');
    end
end
fprintf('\n');

% Get the phase
E2(abs(E2) < 0.001) = 0;
P2 = angle(E2);
p2 = P2(m/2+1,:);
for i = 1:n
    P2(:,i) = unwrap(P2(:,i));
    P2(:,i) = P2(:,i) - P2(m/2+1,i) + p2(i) + pi*x2(:,i).^2./(lambda*Rp(end,i)) + 2*pi*d2(i)/lambda/1000; % Reintroducing the phase from propagation, but scaled down by a factor of 1000 to make the phase curvature more easily seen as the aspect ratio of the simulation is about 1:1000
end
P2 = P2 - P2(m/2+1,round(n/2));

% Interpolate the phase
P = zeros(length(x),n);
for i = 1:n
    P(:,i) = interp1(x2(:,i),P2(:,i),x,'spline',0);
    if mod(i,10) == 0
        fprintf('.');
    end
end
fprintf('\n');


%% Plot the results
% Plot the intensity
figure;
imagesc(1e3*(d2 - fN),1e6*x,I);
title('Intensity');
xlabel('Position relative to focal point [mm]');
ylabel('x [\mum]');
axis equal tight;
set(gca,'YLim',[-1 1]);

% Plot the phase
figure;
imagesc(1e3*(d2 - fN),1e6*x,P,[-1e5 1e5]);
title('Phase');
xlabel('Position relative to focal point [mm]');
ylabel('x [\mum]');

% Plot the phase contour
figure;
contour(1e3*(d2 - fN),1e6*x,P,linspace(-1e5,1e5,41));
title('Phase contour');
xlabel('Position relative to focal point [mm]');
ylabel('x [\mum]');

% Plot the transform order parameter
figure;
plot(1e3*(d2 - fN),sum(a));
title('Tranform order');
xlabel('Position relative to focal point [mm]');
ylabel('Transform order');

% Plot the phase gradient
[Fx,Fy] = gradient(P,mean(diff(d2)),mean(diff(x)));
Fx(I < 1) = 0;
Fy(I < 1) = 0;
k = 10;
figure;
quiver(1e3*(d2(1:k:end)-fN),1e6*x(1:k:end),Fx(1:k:end,1:k:end),Fy(1:k:end,1:k:end));
title('Phase gradient');
xlabel('Position relative to focal point [mm]');
ylabel('x [\mum]');
axis equal tight;
set(gca,'YLim',[-1 1]);


