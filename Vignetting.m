function [sigma_v,r] = Vignetting(R,N,mu,d1,T,f,lambda,yN)
% For use with FrFT and CRLs. Calculate the object plane RMS Gaussian
% vignetting width (sigma_v). The calculation needs the CRL parameters (R)
% radius of curvature of the lenses, (N) number of lenses, (mu) linear
% absorption coefficient, (d1) object-lens distance, (T) lens-lens
% distance, (f) focal length of individual lenses, (lambda) wavelength, and
% (yN) the RMS width of the Gaussian pupil function. The output parameter
% (r) is the residual after fitting, and is only used for debugging.
%
% Example of usage:
% sigma_v = Vignetting(R,N,mu,d1,T,f,lambda,yN)
%
% Debugging, get the residual from fitting:
% [sigma_v,r] = Vignetting(R,N,mu,d1,T,f,lambda,yN)
% 
% Author: Anders F. Pedersen
%

% Set width parameter based on CRL parameters
w0 = sqrt(R/(2*N*mu));

% Set up object parameters
m = 100;
dx = 15*w0;

% Make coordinate
x0 = ((-m/2):(m/2 - 1))';
x0 = x0/m*dx;

% Make object
w = 10*w0;
E0 = rect(x0/w);

% Make the propagation distance and focal length arrays
D = [d1+T/2;T*ones(N-1,1);0];
F = f*ones(N,1);

% Calculate the propagation parameters
R0 = Inf;
s0 = dx/sqrt(m);
[a,~,~,sm,sp,~,~] = FrFT_parameters(D,F,lambda,R0,s0);

% Propagate in N steps with attenuation
att = @(x) exp(-mu*x.^2./(2.*R)); % Using sqrt of attenuation, acts on the amplitude
EN = E0;
xN = x0;
for i = 1:N
    [EN,xN] = propFrFT1(EN,xN,Inf,Inf,sm(i),sp(i),a(i),lambda,0);
    EN = EN.*att(xN);
end

% Fit the CRL exit plane
opt = optimset('TolX',1e-4*w0);
[sigma_v,r] = fminsearch(@(x) fit_fun(x),2*w0,opt);

    function r = fit_fun(sigvig)
        % Apply vignetting
        attvig = exp(-x0.^2./(4.*sigvig.^2)); % Using sqrt of attenuation, acts on the amplitude
        E = E0.*attvig;
        
        % Propagate to CRL exit
        [E,x] = propFrFT1(E,x0,Inf,Inf,sm(1),sp(end-1),sum(a(1:end-1)),lambda,0);
        
        % Apply aperture
        attpsf = exp(-x.^2./(4.*yN.^2));
        E = E.*attpsf;
        
        % Calculate residual
        r = sum(abs(abs(EN).^2 - abs(E).^2).^2);
    end

end

