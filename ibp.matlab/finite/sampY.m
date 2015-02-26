function Ysamp = sampY(X,Y,Z,alpha,epsilon,lambda,p)
pYjointprop = zeros([size(Y) 2]);

N = size(X,2);
K = size(Z,2);

Ysamp = Y;

for(k = 1:size(Ysamp,1))
    for(t = 1:size(Ysamp,2))
        for(a = 0:1)
            oldYkt = Ysamp(k,t);
            Ysamp(k,t) = a;

            Ypriorpart = (p^Ysamp(k,t))*((1-p)^(-Ysamp(k,t)));

            e = Z(:,:)*Ysamp(:,t);
            Xpriorpart = prod(((1-((1-lambda).^e)*(1-epsilon)).^(X(:,t))).*((((1-lambda).^e)*(1-epsilon)).^(1-X(:,t))));

            pykta = Ypriorpart*Xpriorpart;


            Ysamp(k,t) = oldYkt;
            pYjointprop(k,t,a+1) = pykta;
        end
        pY0 = pYjointprop(k,t,1) ./ sum(pYjointprop(k,t,:),3);
        if(pY0 > rand)
            Ysamp(k,t) = 0;
        else
            Ysamp(k,t) = 1;
        end
    end
end
% 
% pY0 = pYjointprop(:,:,1) ./ sum(pYjointprop,3);
% pY1 = 1-pY0;
% % imagesc(Ysamp);
% %
% gi = find(pY0 > rand(size(Y)));
% Ysamp = ones(size(Y));
% Ysamp(gi) = 0;