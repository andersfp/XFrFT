function [E1,x1] = propFrFT1(E0,x0,R0,R1,s0,s1,a,lambda,d)
% Use the fractional Fourier transform (FrFT) to propagate the electric
% field (E0) (1D only) with coordinates (x0) to a detector positioned a
% distance (d) away. The radius of curvature of the initial and final
% planes are given by (R0) and (R1), as well as the scaling factors for the
% initial and final planes (s0) and (s1). The overall transform parameter
% is (a).
% 
% Example of usage:
% [E1,x1] = propFrFT1(E0,x0,R0,R1,s0,s1,a,lambda,d)
%
% The equations have been adopted from:
% Ozaktas, H. M.; Mendlovic, D. Fractional Fourier Optics. J. Opt. Soc. Am.
% A, JOSAA 1995, 12 (4), 743-751.
%
% Author: Anders F. Pedersen
% 

% Make sure the input electric field is a vector
if ~isvector(E0)
    error('Input electric field must be a vector.');
end

% Make sure the input is a column vector
if ~iscolumn(E0)
    E0 = E0.';
end

% Make sure the input coordinate is a vector
if ~isvector(x0)
    error('Input coordinate must be a vector.');
end

% Make sure the input coordinate is a column vector
if ~iscolumn(x0)
    x0 = x0.';
end

% Apply the first chirp (curved phase)
if isinf(R0)
    E1 = E0;
else
    E1 = E0.*exp(-1i.*pi.*x0.^2./(lambda.*R0));
end

% Perform the fractional Fourier transform
E1 = frfft1for(E1,a);

% Calculate the detector plane coordinates
x1 = x0.*(s1./s0);

% Apply the second chirp
if isinf(R1)
    
else
    E1 = E1.*exp(1i.*pi.*x1.^2./(lambda.*R1));
end

% Apply phase factors and amplitude scaling
E1 = E1.*(exp(-1i.*a.*pi./4).*sqrt(s0./s1).*exp(1i.*2.*pi.*d./lambda));


