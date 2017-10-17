function [a,Rm,Rp,sm,sp,gm,gp] = FrFT_parameters(d,f,lambda,Rinit,sinit)
% Calculate the fractional Fourier transform parameters for an arbitrary
% stack of lenses, given the lens separations (d), focal lengths (f), the
% wavelength (lambda), the curvature of the inital plane (Rinit), and the
% scaling parameter of the initial plane (sinit).
%
% Example of usage:
% [a,Rm,Rp,sm,sp,gm,gp] = FrFT_parameters(d,f,lambda,Rinit,sinit)
%
% The equations have been adopted from:
% Ozaktas, H. M.; Mendlovic, D. Fractional Fourier Optics. J. Opt. Soc. Am.
% A, JOSAA 1995, 12 (4), 743-751.
%
% Author: Anders F. Pedersen
% 

% Get number of spacings
n = length(d);

% Prepare arrays
phi = zeros(n,1);
Rm = phi;
Rp = phi;
sm = phi;
sp = phi;
gm = phi;
gp = phi;

% Assign initial values
Rm(1) = Rinit;
sm(1) = sinit;

% Calculate parameters iteratively
for i = 1:n
    gm(i) = 1 + d(i)./Rm(i);
    phi(i) = acot(gm(i).*sm(i).^2./(lambda.*d(i)));
    if d(i) == 0 % Special case for d == 0 to avoid NaN
        sp(i) = sm(i);
        gp(i) = gm(i);
        Rp(i) = Rm(i);
    else
        sp(i) = lambda.*d(i)./sm(i).*csc(phi(i));
        gp(i) = lambda.*d(i)./sp(i).^2.*cot(phi(i));
        Rp(i) = d(i)./(1 - gp(i));
    end
    if sp(i) < 0 % Ensure that all sm and sp are positive. If sp is negative it can be inverted by adding pi to the transform parameter
        sp(i) = -sp(i);
        phi(i) = phi(i) + pi;
    end
    if i < n % Calculate the initial parameters for the next propagation step, unless it is the last step
        sm(i+1) = sp(i);
        Rm(i+1) = 1./(1./Rp(i) - 1./f(i));
    end
end

% Convert the transform degree from radians to order
a = phi./(pi/2);



