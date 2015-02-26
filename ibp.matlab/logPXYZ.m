function lP = logPXYZ(X,Y,Z,alpha,epsilon,lambda,p,Temp)

if(nargin<8)
    Temp =1;
end

nY = prod(size(Y));
nY1 = length(find(Y==1));
nY0 = nY-nY1;

lPY = nY1*log(p)+nY0*log(1-p);

N = size(Z,1);
K = size(Z,2);

HN = sum(1./(1:N));

m_k = sum(Z);
FINITE = 0;
if(FINITE)
    lPZ = sum(gammaln(m_k+alpha/K)+gammaln(N-m_k+1)-gammaln(N+1+alpha/K));
else
    lPZ = K*log(alpha) -alpha*HN + sum(gammaln(m_k)+gammaln(N-m_k+1)-gammaln(N+1));
end

lPX = sum(sum(X.*log(1-(1-lambda).^(Z*Y)*(1-epsilon))+(1-X).*log((1-lambda).^(Z*Y)*(1-epsilon))));

% lP = Temp*(lPX + lPZ+ lPY);
lP = Temp*(lPX + lPZ);

% TLG: we don't want to include lPY, since Y is a big infinite matrix anyway