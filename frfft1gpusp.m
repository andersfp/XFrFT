function res = frfft1gpusp(fc,a)
% Calculate the 1D fractional Fourier transform along the first dimension
% of the input (fc). The transform order is given by the second input (a).
% The input (fc) must have an even number of rows.
% Single precision only. Requires a compatible GPU.
% 
% Example of usage:
% res = frfft1gpusp(fc,a)
%
% The function supports single precision input only.
%
% This implementation is performs the central algorithm on a GPU, giving
% significant speedup over a CPU. The algorithm breaks up the array into a
% size that has the highest FFT performance. This function has been
% optimized for an Nvidia GTX Titan X (Pascal), for which the optimum
% number of elements for FFT is elmax = 2*4.98e7. This number might be
% different other GPUs, so for optimal performance this should be tested
% and the elmax parameter in this function should be changed.
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

% Split the array to optimize FFT computation
%elmax = 2*4.98e7; % This is the number of elements that maximize the FFT performance on an Nvidia Titan X (Pascal) GPU with 12 GB VRAM
elmax = 49152000;
m = floor(elmax/(2^nextpow2(16*N)));
k = ceil(M/m);
i1 = (0:(k - 1))*m + 1;
i2 = (1:k)*m;
i2(end) = M;

% Pre-allocate memory for the result
res = zeros(N,M,'single');

% Make shift array to shift data
shft = single([(4*N/2+1):(4*N) 1:(4*N/2)]);

for i = 1:k
    % Send data to GPU
    fg = gpuArray(fc(:,i1(i):i2(i)));
    
    % Make local version of the a-parameter
    ag = a;
    
    % Interpolate the input function
    fg = bizinter(fg);
    
    % Zeropad the array
    mg = size(fg,2);
    fg = cat(1,zeros(N,mg,'single','gpuArray'),fg,zeros(N,mg,'single','gpuArray'));
    
    % Map a onto the interval 0.5 <= |a| <= 1.5
    if ((ag > 0) && (ag < 0.5)) || ((ag > 1.5) && (ag < 2))
        ag = ag - 1;
        fg(shft,:) = fft(fg(shft,:))/sqrt(4*N);
    elseif ((ag > -0.5) && (ag < 0)) || ((ag > -2) && (ag < -1.5))
        ag = ag + 1;
        fg(shft,:) = ifft(fg(shft,:))*sqrt(4*N);
    end
    
    % Calculate the transform at reduced interval a
    fg = corefrmod2(fg,ag);
    
    % Deinterpolate the result
    fg = fg(N+1:2:3*N,:);
    
    % Return the result to the CPU
    res(:,i1(i):i2(i)) = gather(fg);
end

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
deltax = single(sqrt(N));
phi = a*pi/2;
beta = 1/sin(phi);

% Calculate chirp vectors
x = (-ceil(N/2):fix(N/2)-1).'/deltax;
chrp1 = exp(-1i*pi*tan(phi/2)*x.^2);
t = (-N+1:N-1).'/deltax;
chrp2 = exp(1i*pi*beta*t.^2);

% Get lengths of chirp and fft length
N2 = 2*N - 1;
N3 = 2^nextpow2(N2 + N - 1);

% Zeropad chirp for convolution
chrp2 = cat(1,chrp2,zeros(N3 - N2,1));

% Fourier transform chirp
chrp2 = fft(chrp2);

% Send chirps to the GPU
chrp1 = gpuArray(single(chrp1));
chrp2 = gpuArray(single(chrp2));

% Calculate amplitude
Aphi = single(exp(-1i*(pi*sign(sin(phi))/4-phi/2))/sqrt(abs(sin(phi))));

% Multiply by chirp
fc = bsxfun(@times,fc,chrp1);

% Zeropad array for convolution
fc = cat(1,fc,zeros(N3 - N,M,'single','gpuArray'));

% Perform chirp convolution
fc = fft(fc);
fc = bsxfun(@times,fc,chrp2);
fc = ifft(fc);
fc = fc(N:2*N-1,:);

% Apply amplitude and chirp multiplication
res = bsxfun(@times,fc,chrp1).*(Aphi./deltax);

% Shift array if odd sized array
if rem(N,2) == 1
    res = circshift(res,-1);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xint = bizinter(x)

% Get the number of data points
N = single(size(x,1));
M = single(size(x,2));

% Determine if input is complex, and split real and complex parts
im = 0;
if ~all(isreal(x(:)))
    im = 1;
    imx = imag(x);
    x = real(x);
end;

% Process the real part
xint = bizintercore(x);

% Process the imaginary part
if im == 1
    xmint = bizintercore(imx);
    xint = xint + single(1i)*xmint;
end

% Add core function
    function xint = bizintercore(x2)
        % Add zeros at every other element
        x2 = cat(3,x2,zeros(N,M,'single','gpuArray'));
        x2 = permute(x2,[3 1 2]);
        x2 = reshape(x2,2*N,M);
        
        % Fourier transform the array
        xf = fft(x2);
        
        % Inverse Fourier transform
        if rem(N,2) == 1
            N1 = fix(N/2+1);
            N2 = 2*N - fix(N/2) + 1;
            xf = cat(1,xf(1:N1,:),zeros(N,M,'single','gpuArray'),xf(N2:2*N,:));
            xint = ifft(xf);
            xint = 2*real(xint);
        else
            xf(N/2+1:2*N-N/2,:) = 0;
            xint = ifft(xf);
            xint = 2*real(xint);
        end
    end

end

