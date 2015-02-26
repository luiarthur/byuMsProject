function [Zsamp,Ysamp] = sampZ_finite(X,Y,Z,alpha,epsilon,lambda,p,rows_to_sample)

if(nargin<8)
    rows_to_sample = 1:size(Zsamp,1)

end

Temp = 1;
    
pZjointprop = zeros([size(Z) 2]);
Zsamp = Z;
N = size(Z,1);
K = size(Z,2);

T = size(X,2);
Ysamp = Y;
K = size(Z,2);

zero_inds = find(Z==0);
one_inds = find(Z~=0);
K = size(Zsamp,2);

for(i=rows_to_sample)
    for(k=1:K)
        [Zsamp, Ysamp] = sampleIndividualEntry(X, Zsamp,Ysamp, alpha,epsilon,lambda,p,N, K, T, Temp, i, k);
    end
%     [Zsamp, Ysamp] = sampleNewColumns(X, Zsamp,Ysamp, alpha,epsilon,lambda,p,N, K, T, Temp, i);
%     non_empty_columns = find(sum(Zsamp)>0);
%     Zsamp = Zsamp(:,non_empty_columns);
%     Ysamp = Ysamp(non_empty_columns,:);
%     K = size(Zsamp,2);
end

end

function [Zsamp, Ysamp] = sampleIndividualEntry(X, Zsamp,Ysamp, alpha,epsilon,lambda,p,N, K, T, Temp, i, k)
FINITE = 1;

pZik = zeros(2,1);
mk = sum(Zsamp([1:i-1 i+1:end],k)); %this is m_{-i,k} CORRECT
if (mk > 0 | FINITE)
    oldZik = Zsamp(i,k);
    for(a = 0:1)
        if(FINITE)
            if(a == 1)
                Zpriorpartl = log((mk+alpha/K)/(alpha/K+N));
            else
                Zpriorpartl = log(1-(mk+alpha/K)/(alpha/K+N));
            end
        else
            if(a == 1)
                Zpriorpartl = log(mk/N);
            else
                Zpriorpartl = log(1-(mk/N));
            end
        end


        Zsamp(i,k) = a;

        e = Zsamp(i,1:K)*Ysamp(1:K,:);
        Xpriorpartl = sum(X(i,:).*log(1-((1-lambda).^e)*(1-epsilon)) ...
            +(1-X(i,:)).*log((1-lambda).^e)*(1-epsilon));

        pzikal = Zpriorpartl+Xpriorpartl;


        pZik(a+1) = Temp*pzikal;

    end
    Zsamp(i,k) = oldZik;

    pZ0 = 1./(1+exp(-(pZik(1)-pZik(2))));
    if(pZ0 > rand)
        Zsamp(i,k) = 0;
    else
        Zsamp(i,k) = 1;
    end
else
    %     Zsamp(i,k) = 0;
end
end

function [Zsamp, Ysamp] = sampleNewColumns(X, Zsamp,Ysamp, alpha,epsilon,lambda,p,N, K, T, Temp, i)

% m = sum(Zsamp(setdiff(1:N,i),:));
% inds = find(m~=0);
% Zsamp = Zsamp(:,inds);
% Ysamp = Ysamp(inds,:);
% K = size(Zsamp,2);

m = sum(Zsamp(setdiff(1:N,i),:));
inds = find(m==0);
Zsamp(i,inds) = 0;

newK = 0;
cdfr = rand;
e = Zsamp(i,1:K)*Ysamp(1:K,:);
oneinds = find(X(i,:)==1);
zeroinds = setdiff(1:T,oneinds);
lpnewK = zeros(10,1);
for(newK = 0:9)

    lpXiT = 0;
    lpXiT = sum(log(1-(1-epsilon)*((1-lambda).^e(oneinds))*((1-lambda*p)^newK)));
    lpXiT = lpXiT+sum(log((1-epsilon)*((1-lambda).^e(zeroinds))*((1-lambda*p)^newK)));
    lpK1i = lpXiT -alpha/N + (K+newK)*log(alpha/N) - gammaln(K+newK+1);

    lpnewK(newK+1) = Temp*lpK1i;
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
    Zsampnew([1:size(Zsamp,1)],[1:kplus]) = Zsamp;
    Zsamp = Zsampnew;
    Zsamp(i,kplus+1:end) = 1;

    %         Ysampnew = zeros(size(Ysamp) + [newK 0]);
    Ysampnew = [Ysamp; zeros(newK,T)];

    % TLG: the new values of Y should be drawn jointly from their
    % posterior distribution given Z

    newprobs = zeros(T,newK+1);
    oneinds = find(X(i,:)==1);
    zeroinds = setdiff(1:T,oneinds);
    newprobs = ((1-epsilon)*(1-lambda).^(repmat(e,newK+1,1)+repmat((0:newK)',1,T)));
    newprobs(:,oneinds) = 1 - newprobs(:,oneinds);

    newprobs = newprobs.*repmat((nchoosekv(newK,0:newK).*(p.^(0:newK)).*((1-p).^(newK-(0:newK))))',1,T);
    newprobs = newprobs.^Temp;

    newprobs = newprobs./repmat(sum(newprobs),newK+1,1);
    newprobs = cumsum(newprobs);
%     if(newK>1)
%         disp('fo0')
%     end
    
    for(j=1:T)
        m = min(find(rand<newprobs(:,j)))-1;
        Ysampnew(end-m+1:end,j) = 1;
    end
    
%     m = newprobs < repmat(rand(1,size(newprobs,2)),newK+1,1);
%     Ysampnew((end-newK+1):end,:) = m(1:end-1,:);   % this one is questionable
    %
%             for j = 1:T
%                 newprobs = zeros(1,newK+1);
%                 for m = 0:newK
%                     if (X(i,j) == 1)
%                         newprobs(m+1) = 1-(1-epsilon)*(1-lambda)^(e(j)+m);
%                     else
%                         newprobs(m+1) = (1-epsilon)*(1-lambda)^(e(j)+m);
%                     end
%                     newprobs(m+1) = ((newprobs(m+1)*nchoosek(newK,m)*p^m*(1-p)^(newK-m)));
%                 end
%                 newprobs = newprobs/sum(newprobs);
%                 newprobs = cumsum(newprobs);
%                 m = min(find(rand<newprobs))-1;
%                 Ysampnew(1:m,j) = ones(m,1);
%             end
    Ysamp = Ysampnew;
end
end


function retval = nchoosekv(n,k)
retval = zeros(size(k));
for(i=1:length(k))
    retval(i) = nchoosek(n,k(i));
end
end