% generate binary observations matrix from underlying graphical model then
% run Gibbs, RJMCMC, and SIS samplers to verify that the estimated
% posterior is peaked around the true parameter values

N = 6;                  % number of rows in observation matrix X
T = 100;                % number of columns in observation matrix X
alpha = 2;              % IBP concentration parameter
lambda = .8;            % noisy-or parameters (see paper)
epsilon = .1;           %
p = .4;                 %
num_samples = 20000;    % samples to draw using the Gibbs and RJMCMC sampler
%num_samples = 5000;    % samples to draw using the Gibbs and RJMCMC sampler

% generate test data from the model
[X,Z,Y] = generate_test_data(N,T,alpha,lambda,epsilon,p);

% run the Gibbs sampler (called hyper_sampler because this sampler samples
% the hyperparameters in the model as well
[Z_sample, Y_sample, lP_sample, K_sample, alpha_sample, epsilon_sample, lambda_sample, p_sample] = hyper_sampler(X,num_samples,Y,Z,alpha,epsilon,lambda,p)

% compute some posterior averages
[Ek, EZZt] = inferstats(Z_sample,Z,0);

% run the RJMCMC sampler
[Z_sample, Y_sample, lP_sample, K_sample, alpha_sample, epsilon_sample, lambda_sample, p_sample] = rjmcmc_sampler(X,num_samples,1,Y,Z,alpha,epsilon,lambda,p)

num_particles = 1000;   % number of particles to use in the SIS sampler

% run the SIS sampler
[Zparticles, Yparticles] = particle_filter(X,Z,Y,alpha,epsilon,lambda,p,num_particles)
