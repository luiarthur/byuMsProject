function lP = logPZ_ddibp(Kn,logS,C,alpha,h)
    
    lP = 0;
    for n = 1:length(h)
        lambda = sum(Kn==n);
        lP = lP + log(poisspdf(lambda,alpha/h(n)));
        for k = 1:size(C,2)
            lP = lP + logS(n,C(n,k));
        end
    end