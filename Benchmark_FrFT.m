% Initialization (optional)
clear;
close all;
clc;


%% Benchmark
% Array sizes to test (must be even)
m = logspace(2,3,11).';
m = 2*round(m/2);
n = length(m);

% Number of averages at each size
k = 2;

% Start the parallel pool if not already running
p = gcp('nocreate');
if isempty(p)
    p = parpool;
end

% Transform parameter
a = 0.7;

% Benchmark double precision
T_for_d = zeros(n,k+1);
T_par_d = T_for_d;
T_vec_d = T_for_d;
T_gpu_d = T_for_d;
for i = 1:n
    A = rand(m(i));
    for j = 1:k+1
        tic;
        B = frfft1for(A,a);
        T_for_d(i,j) = toc;
    end
    for j = 1:k+1
        tic;
        C = frfft1par(A,a);
        T_par_d(i,j) = toc;
    end
    for j = 1:k+1
        tic;
        D = frfft1vec(A,a);
        T_vec_d(i,j) = toc;
    end
    for j = 1:k+1
        tic;
        E = frfft1gpu(A,a);
        T_gpu_d(i,j) = toc;
    end
    fprintf('.');
end
fprintf('\n');

% Get average time
t_for_d = mean(T_for_d(:,2:end),2);
t_par_d = mean(T_par_d(:,2:end),2);
t_vec_d = mean(T_vec_d(:,2:end),2);
t_gpu_d = mean(T_gpu_d(:,2:end),2);

% Benchmark single precision
T_for_s = zeros(n,k+1);
T_par_s = T_for_s;
T_vec_s = T_for_s;
T_gpu_s = T_for_s;
for i = 1:n
    A = rand(m(i),'single');
    for j = 1:k+1
        tic;
        B = frfft1for(A,a);
        T_for_s(i,j) = toc;
    end
    for j = 1:k+1
        tic;
        C = frfft1par(A,a);
        T_par_s(i,j) = toc;
    end
    for j = 1:k+1
        tic;
        D = frfft1vec(A,a);
        T_vec_s(i,j) = toc;
    end
    for j = 1:k+1
        tic;
        E = frfft1gpusp(A,a);
        T_gpu_s(i,j) = toc;
    end
    fprintf('.');
end
fprintf('\n');

% Get average time
t_for_s = mean(T_for_s(:,2:end),2);
t_par_s = mean(T_par_s(:,2:end),2);
t_vec_s = mean(T_vec_s(:,2:end),2);
t_gpu_s = mean(T_gpu_s(:,2:end),2);

% Plot the timings
figure;
loglog(m,[t_for_d t_par_d t_vec_d t_gpu_d]);
hold on;
set(gca,'ColorOrderIndex',1);
loglog(m,[t_for_s t_par_s t_vec_s t_gpu_s],'--');
title('Double (single) precision timing');
xlabel('Image size');
ylabel('Runtime [s]');
legend('For','Par','Vec','GPU','Location','NorthWest');

% Save the simulation results
%clear A B C D E;
%save('Benchmark.mat');


