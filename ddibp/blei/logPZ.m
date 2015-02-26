function lP = logPZ(Z,alpha)
    
    Z(:,sum(Z)==0) = [];
    N = size(Z,1);
    K = size(Z,2);
    
    HN = sum(1./(1:N));
    
    m_k = sum(Z);
    
    lP = K*log(alpha) -alpha*HN + sum(gammaln(m_k)+gammaln(N-m_k+1)-gammaln(N+1));