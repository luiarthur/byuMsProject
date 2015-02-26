function Zsamp = sampZ(X,Z,sigma_x, sigma_A, alpha)


N = size(X,1);
K = size(Z,2);
D = size(X,2);

Zsamp = Z;

for i=1:N
    for k=1:K
        Zsamp = sampleIndividualEntry(X, Zsamp, sigma_x, sigma_A,N, K, D, i, k);
    end
    Zsamp = sampleNewColumns(X, Zsamp, sigma_x, sigma_A, alpha ,N, K, D, i);
    Zsamp = Zsamp(:,sum(Zsamp)>0);
    K = size(Zsamp,2);
end

end

function Zsamp = sampleIndividualEntry(X, Zsamp,sigma_x, sigma_A, N, K, D, i, k)

pZik = zeros(2,1);
mk = sum(Zsamp([1:i-1 i+1:end],k)); %this is m_{-i,k} CORRECT

if (mk > 0)
    oldZik = Zsamp(i,k);
    for a = 0:1

        if(a == 1)
            logZprior = log(mk/N);
        else
            logZprior = log(1-(mk/N));
        end

        Zsamp(i,k) = a;

        logZpart = (D/2)* log (det(Zsamp'*Zsamp + (sigma_x^2)/(sigma_A^2)*eye(K)));
        logExp = -1/(2*sigma_x^2) * trace(X'*(eye(N)- Zsamp*inv(Zsamp'*Zsamp + (sigma_x^2)/(sigma_A^2)*eye(K))*Zsamp')*X);


        pzikal = logExp-logZpart+logZprior;


        pZik(a+1) = pzikal;

    end
    Zsamp(i,k) = oldZik;

    pZ0 = 1./(1+exp(-(pZik(1)-pZik(2))));
    if(pZ0 > rand)
        Zsamp(i,k) = 0;
    else
        Zsamp(i,k) = 1;
    end
else
    %     do nothing
end
end

function Zsamp = sampleNewColumns(X, Zsamp,sigma_x, sigma_A, alpha,N, K, D, i)

m = sum(Zsamp(setdiff(1:N,i),:));
Zsamp(i,m==0) = 0;

cdfr = rand;

lpnewK = zeros(10,1);

for newK = 0:9

    logZ = (N*D/2) * log(2*pi) + ((N-K)*D) * log(sigma_x) +  (K * D) * log(sigma_A) * (D/2)* log (det(Zsamp'*Zsamp + (sigma_x^2)/(sigma_A^2)*eye(K)));
    logExp = -1/(2*sigma_x^2) * trace(X'*(eye(N)- Zsamp*inv(Zsamp'*Zsamp + (sigma_x^2)/(sigma_A^2)*eye(K))*Zsamp')*X);
    logLike = logExp - logZ;
    
    lpK1i = logLike -alpha/N + (K+newK)*log(alpha/N) - gammaln(K+newK+1);

    lpnewK(newK+1) = lpK1i;
end
logmax = max(lpnewK);

pdf = exp(lpnewK-logmax);
pdf = pdf/sum(pdf);
cdf = pdf(1);
newK = 0;
ii=1;
while(cdf<cdfr)
    ii=ii+1;
    cdf=cdf+pdf(ii);
    newK = newK+1;
end

if(newK>0)
    kplus = size(Zsamp,2);
    Zsampnew = zeros(size(Zsamp) +[0 newK]);
    Zsampnew(1:size(Zsamp,1),1:kplus) = Zsamp;
    Zsamp = Zsampnew;
    Zsamp(i,kplus+1:end) = 1;
end

end
