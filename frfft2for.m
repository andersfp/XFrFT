function F2D = frfft2for(f2D,a)
% Perform a 2D fractional Fourier transform of the input (f2D). The 
% transform order (a) should be specified as [ax ay]. It is based on the 1D
% version frfft1for.
% 
% Example of usage:
% F2D = frfft2for(f2D,a)
% 
% Author: Anders F. Pedersen
% 

% Get the individual transform orders
ax = a(1);
ay = a(2);

% Perform the transform in the y-direction (columns, vertical)
F2D = frfft1for(f2D,ay);

% Perform the transform in the x-direction (rows, horizontal)
F2D = F2D.';
F2D = frfft1for(F2D,ax);
F2D = F2D.';

