function out = rect(x)
% Rectangle function. Returns 1 for |x| <= 1/2, 0 otherwise.
% 
% Example of usage:
% out = rect(x)
%

out = double(abs(x) <= 1/2);

