% establish 

addpath ..
addpath .
addpath ../../../utilities

load FACEDATA.mat

% [Z_sample, lP_sample, K_sample, alpha_sample, sigma_x_sample, sigma_A_sample] =hyper_sampler(X,100,Z_plus,2,sigma_x, 1)
% 
% Z_last = Z_sample{end};Aest = inv(Z_last'*Z_last)*Z_last'*X
% 
% figure
% num_figs = ceil(sqrt(size(Aest,1)));
% for i=1:size(Aest,1)
%     subplot(num_figs,num_figs,i)
%     imagesc(reshape(Aest(i,:),6,6))
% end

Zparticles = particle_filter(zmX,10,4,.25, .01)

Z_last_pf = Zparticles{1};Aestpf = pinv(Z_last_pf'*Z_last_pf)*Z_last_pf'*X(1:n,:)
figure
num_figs = ceil(sqrt(size(Aestpf,1)));
for i=1:size(Aestpf,1)
    subplot(num_figs,num_figs,i)
    imagesc(reshape(Aestpf(i,:),19,19))
end

save FACEDATA-alpha-10-100-particles.mat