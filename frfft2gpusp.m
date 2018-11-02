function F2D = frfft2gpusp(f2D,a)
% Perform a 2D fractional Fourier transform of the input (f2D). The 
% transform order (a) should be specified as [ax ay]. It is based on the 1D
% version frfft1gpusp. The first dimension is treated as y and the second
% dimension as x, any further dimensions will not be transformed, so e.g. a
% 3D array will be treated as a stack of 2D images.
% Single precision only. Requires a compatible GPU.
% 
% Example of usage:
% F2D = frfft2gpusp(f2D,a)
% 
% Author: Anders F. Pedersen
% 

% Get the individual transform orders
ax = a(1);
ay = a(2);

% Convert input to single precision
f2D = single(f2D);

% Perform the transform in the y-direction (columns, vertical)
F2D = frfft1gpusp(f2D,ay);

% Perform the transform in the x-direction (rows, horizontal)
n = ndims(F2D);
ii = 1:n;
ii(1:2) = ii(2:-1:1);
F2D = permute(F2D,ii);
F2D = frfft1gpusp(F2D,ax);
F2D = permute(F2D,ii);

