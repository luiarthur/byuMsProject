function Ysamp = sampY_newrows_only(X,Y,Z,alpha,epsilon,lambda,p,start_row)

if(nargin<8)
    Temp =1 ;
end

pYjointprop = zeros([size(Y) 2]);

N = size(X,2);
K = size(Z,2);

Ysamp = Y;

for(t = 1:size(Ysamp,2))
    for(k = start_row:size(Ysamp,1))
        for(a = 0:1)
            oldYkt = Ysamp(k,t);
            Ysamp(k,t) = a;

            Ypriorpartl = Ysamp(k,t)*log(p)+(1-Ysamp(k,t))*log(1-p);

            e = Z(:,:)*Ysamp(:,t);
            Xpriorpartl = sum(X(:,t).*log(1-((1-lambda).^e)*(1-epsilon)) ...
                + (1-X(:,t)).*log(((1-lambda).^e)*(1-epsilon))); %remember that we need p(x) not p(x==1)

            pyktal = Ypriorpartl+Xpriorpartl;


            Ysamp(k,t) = oldYkt;
            pYjointprop(k,t,a+1) = pyktal; %temp applied here
        end
        pY0 = 1/(1+exp(-(pYjointprop(k,t,1)-pYjointprop(k,t,2))));  % <-- good trick 1/(1+e^(lp1/lp0) = p0/(p0+p1)
        if(pY0 > rand)
            Ysamp(k,t) = 0;
        else
            Ysamp(k,t) = 1;
        end
    end
end
