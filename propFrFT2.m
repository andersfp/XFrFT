function [E1,X1,Y1] = propFrFT2(E0,X0,Y0,R0,R1,s0,s1,a,lambda,d,method)
% Use the fractional Fourier transform (FrFT) to propagate the electric
% field (E0) (2D or ND stacks of 2D planes) with coordinates (X0 and Y0) to
% a detector positioned a distance (d) away. The radius of curvature of the
% initial and final planes are given by (R0) and (R1), as well as the 
% scaling factors for the initial and final planes (s0) and (s1). The
% overall transform parameter is (a). The parameters R0, R1, s0, s1, and a
% maybe be given as scalars if identical for the x- and y-directions, and
% may be given as [px py] is they are different. The last input (method) is
% used to decide which algorithm to be used, either 'gpu', 'gpusp', 'vec',
% 'par', or 'for'.
% 
% Example of usage:
% [E1,X1,Y1] = propFrFT2(E0,X0,Y0,R0,R1,s0,s1,a,lambda,d,method)
%
% The equations have been adopted from:
% Ozaktas, H. M.; Mendlovic, D. Fractional Fourier Optics. J. Opt. Soc. Am.
% A, JOSAA 1995, 12 (4), 743-751.
%
% Author: Anders F. Pedersen
% 

% Assign default method if not specified
if nargin == 10
    method = 'vec';
end

% Make sure that x0 is a row vector or matrix
if isvector(X0)
    if ~isrow(X0)
        X0 = X0.';
    end
elseif ~ismatrix(X0)
    error('x0 coordinate must be a vector or matrix.');
end

% Make sure that y0 is a column vector or matrix
if isvector(Y0)
    if ~iscolumn(Y0)
        Y0 = Y0.';
    end
elseif ~ismatrix(Y0)
    error('y0 coordinate must be a vector or matrix.');
end

% Assign radius of curvature for origin
if isscalar(R0)
    R0x = R0;
    R0y = R0;
elseif isvector(R0) && length(R0) == 2
    R0x = R0(1);
    R0y = R0(2);
else
    error('R0 must have either 1 or 2 scalars.');
end

% Assign radius of curvature for detector
if isscalar(R1)
    R1x = R1;
    R1y = R1;
elseif isvector(R1) && length(R1) == 2
    R1x = R1(1);
    R1y = R1(2);
else
    error('R1 must have either 1 or 2 scalars.');
end

% Assign scaling parameter for origin
if isscalar(s0)
    s0x = s0;
    s0y = s0;
elseif isvector(s0) && length(s0) == 2
    s0x = s0(1);
    s0y = s0(2);
else
    error('s0 must have either 1 or 2 scalars.');
end

% Assign scaling parameter for detector
if isscalar(s1)
    s1x = s1;
    s1y = s1;
elseif isvector(s1) && length(s1) == 2
    s1x = s1(1);
    s1y = s1(2);
else
    error('s1 must have either 1 or 2 scalars.');
end

% Assign transform parameter for origin
if isscalar(a)
    ax = a;
    ay = a;
elseif isvector(a) && length(a) == 2
    ax = a(1);
    ay = a(2);
else
    error('a must have either 1 or 2 scalars.');
end

% Calculate the detector plane coordinates
X1 = X0.*(s1x./s0x);
Y1 = Y0.*(s1y./s0y);

% Apply the first chirp (curved phase)
if isinf(R0x) && isinf(R0y)
    E1 = E0;
elseif isinf(R0x)
    E1 = E0.*exp(-1i.*pi.*Y0.^2./(lambda.*R0y));
elseif isinf(R0y)
    E1 = E0.*exp(-1i.*pi.*X0.^2./(lambda.*R0x));
else
    E1 = E0.*exp(-1i.*pi.*X0.^2./(lambda.*R0x)).*exp(-1i.*pi.*Y0.^2./(lambda.*R0y));
end

% Perform the fractional Fourier transform
switch method
    case 'gpu'
        E1 = frfft2gpu(E1,[ax ay]);
    case 'gpusp'
        E1 = frfft2gpusp(E1,[ax ay]);
    case 'par'
        E1 = frfft2par(E1,[ax ay]);
    case 'for'
        E1 = frfft2for(E1,[ax ay]);
    case 'vec'
        E1 = frfft2vec(E1,[ax ay]);
    otherwise
        error('The method must be either ''gpu'', ''gpusp'', ''vec'', ''par'', or ''for''.');
end

% Apply the second chirp
if isinf(R1x) && isinf(R1y)
    
elseif isinf(R1x)
    E1 = E1.*exp(1i.*pi.*Y1.^2./(lambda.*R1y));
elseif isinf(R1y)
    E1 = E1.*exp(1i.*pi.*X1.^2./(lambda.*R1x));
else
    E1 = E1.*exp(1i.*pi.*X1.^2./(lambda.*R1x)).*exp(1i.*pi.*Y1.^2./(lambda.*R1y));
end

% Apply phase factors and amplitude scaling
E1 = E1.*(exp(-1i.*ax.*pi./4).*exp(-1i.*ay.*pi./4).*sqrt(s0x./s1x).*sqrt(s0y./s1y).*exp(1i.*2.*pi.*d./lambda));


