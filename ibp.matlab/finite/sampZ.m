function Zsamp = sampZ(X,Y,Z,alpha,epsilon,lambda,p)

pZjointprop = zeros([size(Z) 2]);
Zsamp = Z;
N = size(Z,1);
K = size(Z,2);


for(k = 1:size(Zsamp,2))
    for(i = 1:size(Zsamp,1))
        for(a = 0:1)


            if(a == 1)
                Zpriorpart = (sum(Z([1:i-1 i+1:end],k))+alpha/K)/(alpha/K+N);
            else
                Zpriorpart = 1-(sum(Z([1:i-1 i+1:end],k))+alpha/K)/(alpha/K+N);
            end

            oldZik = Zsamp(i,k);
            Zsamp(i,k) = a;
%Zpriorpart = 1;
%             mk = sum(Zsamp(:,k));


            e = Zsamp(i,:)*Y(:,:);
            Xpriorpart = prod(((1-((1-lambda).^e)*(1-epsilon)).^X(i,:)).*(((1-lambda).^e)*(1-epsilon)).^(1-X(i,:)));

            pzika = Zpriorpart*Xpriorpart;
            if(pzika==0)
                warning('pzika == 0 problem')
            end


            Zsamp(i,k) = oldZik;
            pZjointprop(i,k,a+1) = pzika;

        end
        pZ0 = pZjointprop(i,k,1) ./ sum(pZjointprop(i,k,:),3);
        if(pZ0 > rand)
            Zsamp(i,k) = 0;
        else
            Zsamp(i,k) = 1;
        end
    end
end
% 
% pZ0 = pZjointprop(:,:,1) ./ sum(pZjointprop,3);
% pZ1 = 1-pZ0;
% imagesc(Zsamp);
%
% gi = find(pZ0 > rand(size(Z)));
% Zsamp = ones(size(Z));
% Zsamp(gi) = 0;