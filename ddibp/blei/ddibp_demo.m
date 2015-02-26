clear

N = 10;
M = 12;
K = 8;
nIter = 100;

sx = 1;
sw = 1;
alpha = 10;

D = zeros(N);
f = @(d) ones(size(d));

W = randn(K,M);
Z = ones(N,K);
X = Z*W + randn(N,M)*0.01;

samples_ddibp = ddibp_mcmc(X,D,f,nIter,sx,sw,alpha);
Z0 = samples_ddibp(1).Z;
samples_ibp = ibp_mcmc(X,nIter,sx,sw,alpha,[],Z0);

score(:,1) = [samples_ddibp.score];
score(:,2) = [samples_ibp.score];

figure;
subplot(2,2,1); imagesc(cannonize(samples_ddibp(end).Z)); title('DD-IBP','FontSize',15,'FontWeight','Bold');
xlabel('Feature','FontSize',15); ylabel('Object','FontSize',15)
subplot(2,2,2); imagesc(cannonize(samples_ibp(end).Z)); title('IBP','FontSize',15,'FontWeight','Bold')
xlabel('Feature','FontSize',15); ylabel('Object','FontSize',15)
subplot(2,2,[3 4]); plot(score,'LineWidth',1.5); legend({'DD-IBP' 'IBP'},'Location','SouthEast','FontSize',15)
ylabel('Log-likelihood','FontSize',15); xlabel('Iteration','FontSize',15);