function lP = logPXYZ(X,Y,Z,alpha,epsilon,lambda,p)

nY = prod(size(Y));
nY1 = length(find(Y==1));
nY0 = nY-nY1;

lPY = nY1*log(p)+nY0*log(1-p);

N = size(Z,1);
K = size(Z,2);

m_k = sum(Z);

lPZ = sum(gammaln(m_k+alpha/K)+gammaln(N-m_k+1)-gammaln(N+1+alpha/K));

lPX = sum(sum(X.*log(1-(1-lambda).^(Z*Y)*(1-epsilon))+(1-X).*log((1-lambda).^(Z*Y)*(1-epsilon))));

lP = lPX + lPZ + lPY;