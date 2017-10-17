function [sigma_D,sigma_a,sigma_v,gamma,yN] = CRL_Parameters_2(N,R,mu,f,phi,d1)
% Calculate the RMS width of the effective aperture (sigma_D), the RMS
% angular acceptance (sigma_a), the RMS vignetting (sigma_v), the angular
% displacement parameter (gamma), and the ray location at the exit of the
% lens stack (yN) of a compund refractive lens (CRL). As input the number
% of lenses (N), radius of curvature (R), linear absorption coefficient
% (mu), individual lens focal length (f), optical focusing power (phi), and
% the sample-lens distance (d1) are needed.
% 
% Example of usage:
% [sigma_D,sigma_a,sigma_v,gamma,yN] = CRL_Parameters_2(N,R,mu,f,phi,d1)
% 
% The formulas are from: Simons, H.; Ahl, S. R.; Poulsen, H. F.;
% Detlefs, C. Simulating and Optimizing Compound Refractive Lens-Based
% X-Ray Microscopes. J Synchrotron Rad, J Synchrotron Radiat 2017, 24 (2), 
% 392-401.
% 
% Author: Anders F. Pedersen
% 

% Supporting parameters
S1 = N.*sinc(2.*N.*phi./pi);
S2 = sin(N.*phi).^2./phi;
D = d1./(f.*phi);

% sigma_D (aperture)
sigma_D = sqrt(R./mu)./sqrt(N + S1);

% sigma_a (angular acceptance)
sigma_a = sqrt(R./mu)./(f.*phi)./sqrt(N.*(1 + D.^2) - (1 - D.^2).*S1 + 2.*D.*S2);

% sigma_v (vignetting)
sigma_v = sqrt(R./mu).*sqrt((N.*(1 + D.^2) - (1 - D.^2).*S1 + 2.*D.*S2)./(N.^2 - S1.^2 - S2.^2));

% gamma (angular conversion)
gamma = (1./(f.*phi)).*(D.*(N + S1) + S2)./(N.*(1 + D.^2) - (1 - D.^2).*S1 + 2.*D.*S2);

% y_N (ray location at CRL exit)
yN = d1.*sigma_a.*(cos(N.*phi) + phi./2.*sin(N.*phi)) + f.*phi.*sigma_a.*(sin(N.*phi) - phi./2.*cos(N.*phi));

