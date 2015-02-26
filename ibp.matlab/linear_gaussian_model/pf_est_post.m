function pf_est_post(fn, zm, sigma_X, sigma_A, alpha)
% establish

disp(['Loading file name = ' fn])
if(ischar(zm))
zm = str2num(zm)
end

if(nargin >2)
    if(ischar(sigma_X))
        sigma_X = str2num(sigma_X);
    end
        if(ischar(sigma_A))
        sigma_A = str2num(sigma_A);
        end
        if(ischar(alpha))
        alpha = str2num(alpha);
    end
else
    sigma_X = .25;
    sigma_A = .01;
    alpha = 4;
end

addpath ..
addpath .
addpath ../../../utilities

load( fn)


if(zm == 1)
    disp(['Zeroing mean'])
    zmX = X-repmat(mean(X),size(X,1),1);
else
    disp(['Leaving data alone'])
    zmX = X;
end

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

Zparticles = particle_filter_for_faces(zmX,100,alpha,sigma_X, sigma_A, ['./tests/pf-big/' fn '_' num2str(zm) '_100_']);

fn = ['./tests/pf-big/' fn '-' num2str(zm) '_100_' num2str(n)];
save(fn)

Z_last_pf = Zparticles{1};
Aestpf = pinv(Z_last_pf'*Z_last_pf)*Z_last_pf'*zmX;
im_size = sqrt(size(zmX,2));
figure
num_figs = ceil(sqrt(size(Aestpf,1)));
for i=1:size(Aestpf,1)
    subplot(num_figs,num_figs,i)
    imagesc(reshape(Aestpf(i,:),im_size,im_size))
end
colormap gray

print(gcf, '-depsc2',[fn '.eps'])

fn = ['./tests/pf-big/' fn '-' num2str(zm) '_100_' num2str(n)];
save(fn)

