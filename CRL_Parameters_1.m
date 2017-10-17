function [f,phi,fN] = CRL_Parameters_1(R,T,N,delta)
% Calculate the individual focal length (f), the optical focusing power
% (phi), and the effective focal length (fN) of a compound refractive lens 
% (CRL), given the physical parameters radius (R), distance between lens
% centers (T), number of lenses (N), and the refractive index decrement
% (delta).
% 
% Example of usage:
% [f,phi,fN] = CRL_Parameters_1(R,T,N,delta)
% 
% The formulas are from: Simons, H.; Ahl, S. R.; Poulsen, H. F.;
% Detlefs, C. Simulating and Optimizing Compound Refractive Lens-Based
% X-Ray Microscopes. J Synchrotron Rad, J Synchrotron Radiat 2017, 24 (2), 
% 392-401.
% 
% Author: Anders F. Pedersen
% 

% Single lens focal length
f = R./(2*delta);

% Optical power per length
phi = sqrt(T./f);

% Effective focal length
fN = f.*phi.*cot(N.*phi);

