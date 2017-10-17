function [delta,mu] = Be_Prop(E)
% Return the delta and mu of Be for each energy in E. The function
% interpolates the data in "Be_delta.txt" and "Be_mu.txt", which has data
% in the energy range 1 keV < E < 30 keV.
% 
% Example of usage:
% [delta,mu] = Be_Prop(E)
% 
% The data is from: http://henke.lbl.gov/optical_constants/
% 
% Author: Author: Anders F. Pedersen
% 

% Load the delta parameter
dat = dlmread('Be_delta.txt');
Ed = dat(:,1);
del = dat(:,2);
%bet = dat(:,3);

% Load the mu parameter
dat = dlmread('Be_mu.txt');
Em = dat(:,1);
m = dat(:,2);

% Interpolate delta
delta = interp1(Ed,del,E,'linear');

% Interpolate mu
mu = interp1(Em,m,E,'linear');

% Convert mu to inverse meters
mu = 1e6./mu;

% Generate warning if energy outside range
if any(isnan(delta))
    warning('Energy outside know range of delta.');
end
if any(isnan(mu))
    warning('Energy outside know range of mu.');
end

