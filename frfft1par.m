function res = frfft1par(fc,a)
% Calculate the 1D fractional Fourier transform along the first dimension
% of the input (fc). The transform order is given by the second input (a).
% The input (fc) must have an even number of rows.
% Requires Parallel Toolbox.
% 
% Example of usage:
% res = frfft1par(fc,a)
%
% The function supports double and single precision inputs.
%
% This implementation is based on a parallel for-loop over each column of 
% the input, and therefore require the Parallel Toolbox.
%
% The basic algorithm is based on an implementation by M. A. Kutay, based
% on the following works:
% Haldun M. Ozaktas, Orhan Arikan, M. Alper Kutay, and Gozde Bozdagi,
% Digital computation of the fractional Fourier transform,
% IEEE Transactions on Signal Processing, 44:2141--2150, 1996.
% Haldun M. Ozaktas, Zeev Zalevsky, and M. Alper Kutay,
% The Fractional Fourier Transform with Applications in Optics and
% Signal Processing, Wiley, 2000, chapter 6, page 298.
%
% Some suggestions A. Bultheel and H. E. M. Sulbaran have been used:
% Bultheel, A.; Martinez Sulbaran, H. E. Computation of the Fractional 
% Fourier Transform. Applied and Computational Harmonic Analysis 2004, 
% 16 (3), 182-202.
%
% Significant speedups and adaptation to 2D array have been made by Anders
% F. Pedersen.
% 
% Author: Anders F. Pedersen
%

% Number of data points in the transform direction
N = size(fc,1);

% Check that the input length is even
if mod(N,2) == 1
    error('Length of the input vector should be even.');
end

% Change a to the interval [-2:2[
a = mod(a + 2,4) - 2;

% Deal with special cases
if a == 0
    res = fc;
    return
elseif a == 2 || a == -2
    res = flip(fc,1);
    return
end

% Reshape ND array to 2D
s = size(fc);
fc = reshape(fc,s(1),prod(s(2:end)));

% Number of data points in the non-transform direction
M = size(fc,2);

% Interpolate the input function
fc = bizinter(fc);
fc = cat(1,zeros(N,M,'like',fc),fc,zeros(N,M,'like',fc));

% Map a onto the interval 0.5 <= |a| <= 1.5
if ((a > 0) && (a < 0.5)) || ((a > 1.5) && (a < 2))
    a = a - 1;
    res = ifftshift(fft(fftshift(fc,1)),1)/sqrt(4*N);
elseif ((a > -0.5) && (a < 0)) || ((a > -2) && (a < -1.5))
    a = a + 1;
    res = ifftshift(ifft(fftshift(fc,1)),1)*sqrt(4*N);
else
    res = fc;
end

% Calculate the transform at reduced interval a
res	= corefrmod2(res,a);

% Deinterpolate the result
res = res(N+1:2:3*N,:);

% Transform output from 2D to ND
res = reshape(res,s);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = corefrmod2(fc,a)
% Core function for computing the fractional Fourier transform.
% Valid only when 0.5 <= abs(a) <= 1.5
% Decomposition used:
%   chirp mutiplication - chirp convolution - chirp mutiplication

% Calculate scalar parameters
N = size(fc,1);
M = size(fc,2);
deltax = sqrt(N);
phi = a*pi/2;
beta = 1/sin(phi);

% Calculate chirp vectors
x = (-ceil(N/2):fix(N/2)-1).'/deltax;
chrp1 = exp(-1i*pi*tan(phi/2)*x.^2);
t = (-N+1:N-1).'/deltax;
chrp2 = exp(1i*pi*beta*t.^2);
clear x;
clear t;
chrp1 = cast(chrp1,'like',fc);
chrp2 = cast(chrp2,'like',fc);

% Get lengths of chirp and fft length
N2 = 2*N - 1;
N3 = 2^nextpow2(N2 + N - 1);

% Zeropad chirp for convolution
chrp2 = cat(1,chrp2,zeros(N3 - N2,1,'like',chrp2));

% Fourier transform chirp
chrp2 = fft(chrp2);

% Calculate amplitude
Aphi = exp(-1i*(pi*sign(sin(phi))/4-phi/2))/sqrt(abs(sin(phi)));

% Run the central multiply-colnvolute-multiply algorithm
res = zeros(N,M,'like',fc);
parfor i = 1:M
    % Multiply by chirp
    fci = fc(:,i).*chrp1;
    
    % Zeropad array for convolution
    fci = cat(1,fci,zeros(N3 - N,1,'like',fci));
    
    % Perform chirp convolution
    fci = ifft(fft(fci).*chrp2);
    fci = fci(N:2*N-1,:);
    
    % Apply amplitude and chirp multiplication
    res(:,i) = (chrp1.*fci).*(Aphi./deltax);
end

% Shift array if odd sized array
if rem(N,2) == 1
    res = circshift(res,-1);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xint = bizinter(x)

% Get the number of data points
N = size(x,1);
M = size(x,2);

% Determine if input is complex, and split real and complex parts
im = 0;
if ~all(isreal(x(:)))
    im = 1;
    imx = imag(x);
    x  = real(x);
end;

% Process the real part
xint = bizintercore(x);

% Process the imaginary part
if im == 1
    xmint = bizintercore(imx);
    xint = xint + 1i*xmint;
end

% Add core function
    function xint = bizintercore(x2)
        xint = zeros(2*N,M,'like',x2);
        if rem(N,2) == 1
            N1 = fix(N/2+1);
            N2 = 2*N - fix(N/2) + 1;
            parfor i = 1:M
                xt = cat(1,x2(:,i).',zeros(1,N,'like',x2(:,i)));
                xt = xt(:);
                xf = fft(xt);
                xf = cat(1,xf(1:N1,:),zeros(N,1,'like',xf),xf(N2:2*N,:));
                xint(:,i) = 2*real(ifft(xf));
            end
        else
            parfor i = 1:M
                xt = cat(1,x2(:,i).',zeros(1,N,'like',x2(:,i)));
                xt = xt(:);
                xf = fft(xt);
                xf(N/2+1:2*N-N/2) = 0;
                xint(:,i) = 2*real(ifft(xf));
            end
        end
    end
end

