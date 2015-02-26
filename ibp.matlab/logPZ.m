function lP = logPZ(Z,alpha,epsilon,lambda,p,FINITE)

if(nargin < 6)
    FINITE = 0;
end

N = size(Z,1);
K = size(Z,2);

HN = sum(1./(1:N));

m_k = sum(Z);

if(FINITE)
    lP = K*log(alpha/K)+sum(gammaln(m_k+alpha/K)+gammaln(N-m_k+1)-gammaln(N+1+alpha/K));
else
    lP = K*log(alpha) -alpha*HN + sum(gammaln(m_k)+gammaln(N-m_k+1)-gammaln(N+1));
end