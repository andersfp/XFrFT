% Initialization (optional)
clear;
close all;
clc;


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


%% Define the object
% Number of pixels on each side (must be even)
m = [10000 1000 500 200 100];

% Field-of-view
dx = 100e-6;

% Plotting markers
pm = '-ox+*sdv^<>ph';

% Prepare figure
figure;
set(gcf,'Position',[300 250 1400 600]);

% Run simulation
n = length(m);
for i = 1:n
% Make coordinates
x0 = ((-m(i)/2):(m(i)/2 - 1)).'/m(i)*dx;

% Make object
E0 = sin(x0*2*pi/5e-6).^2.*exp(-(x0 + 15e-6).^2./(2*(10e-6).^2)) + cos(x0*2*pi/5e-6).^2.*exp(-(x0 - 15e-6).^2./(2*(10e-6).^2));

% Set the object plane curvature and scaling parameter
R0 = Inf;
s0 = dx/sqrt(m(i));

% Calculate the propagation parameters
[a,Rm,Rp,sm,sp,gm,gp] = FrFT_parameters(D,F,lambda,R0,s0);

% Check that the product of the g-parameters is in the range 0<=gm*gp<=1
if any(0 > gm.*gp) || any(gm.*gp > 1)
    warning('Product gm*gp is outside limits.');
end

% Propagate in 2 steps
E1 = E0.*sqrt(exp(-x0.^2./(2*sigma_v.^2))); % Apply the vignetting (sqrt(...) due to amplitude)
[E1,x1] = propFrFT1(E1,x0,Inf,Inf,sm(1),sp(end-1),sum(a(1:end-1)),lambda,sum(D(1:end-1))); % Propagate to the end of the CRL in 1 transform by summing the a-parameter
E1 = E1.*sqrt(exp(-x1.^2./(2*sigma_p.^2))); % Apply the effective pupil (sqrt(...) due to amplitude)
[E2,x2] = propFrFT1(E1,x1,Inf,Inf,sm(end),sp(end),a(end),lambda,D(end)); % Propagate to the detector plane

% Get the intensities
I2 = abs(E2).^2;

% Plot the detector plane with effective vignetting and pupil
plot(1e6*x2,I2,pm(i),'LineWidth',1.3);
hold on;
end

% Set plot options
set(gca,'XLim',[-500 500],'FontSize',14,'OuterPosition',[0 0 1 1],'XTick',-500:100:500,'box','on');
xlabel('x [\mum]');
ylabel('Intensity [a.u.]');


% Make a legend
p = 1e6*dx./m;
p = cellfun(@(x) [num2str(x,3) ' \mum'],num2cell(p),'UniformOutput',false);
legend(p);

% Save the figure
s = 0;
if s == 1
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'Position',[pos(1) pos(2) pos(3) pos(4)]);
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print('Figure_S2.pdf','-dpdf','-r0');
    print('Figure_S2.png','-dpng','-r600');
end


