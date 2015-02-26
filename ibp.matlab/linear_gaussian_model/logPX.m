function lP = logPX(X,Z,sigma_x, sigma_A, alpha)

Z_plus = Z;
K_plus = size(Z,2);
N = size(X,1);
D = size(X,2);
% 
% m_k = sum(Z);
% HN = sum(1./(1:N));
% logZPrior = K_plus*log(alpha) -alpha*HN + sum(gammaln(m_k)+gammaln(N-m_k+1)-gammaln(N+1));

% eqn. 58 of tutorial
% Z = (2*pi)^(N*D/2) * sigma_x^((N-K_plus)*D) * sigma_A^(K_plus * D) * det(Z_plus'*Z_plus + (sigma_x^2)/(sigma_A^2)*eye(K_plus))^(D/2);
% exp(-1/(2*sigma_x^2) * trace(X'*(eye(N)- Z_plus*inv(Z_plus'*Z_plus +
% (sigma_x^2)/(sigma_A^2)*eye(K_plus))*Z_plus')*X))


logZ = (N*D/2) * log(2*pi) + ((N-K_plus)*D) * log(sigma_x) +  (K_plus * D) * log(sigma_A) + (D/2)* log (det(Z_plus'*Z_plus + (sigma_x^2)/(sigma_A^2)*eye(K_plus)));
logExp = -1/(2*sigma_x^2) * trace(X'*(eye(N)- Z_plus*inv(Z_plus'*Z_plus + (sigma_x^2)/(sigma_A^2)*eye(K_plus))*Z_plus')*X);

logLike = logExp - logZ;

lP = logLike ;%+ logZPrior;
% if(lP>0)
%     keyboard;
% end