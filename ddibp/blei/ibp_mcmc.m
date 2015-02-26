function samples = ibp_mcmc(X,nIter,sx,sw,alpha,savefile,Z)
    
    %############ INITIALIZATION #############
    
    %--------- parameters --------------%
    N = size(X,1);
    Hn = sum(1./(1:N));
    na = isnan(X);
    missing = any(na(:));
    if missing
        X(na) = randn*sx;
    end
    
    %---------- initialize results structure -----------%
    disp(num2str(1))
    if nargin < 7
        Z = [ones(N,10) zeros(N,1000)];
        samples.Z = Z;
    end
    samples.score = lik(X,Z,sx,sw) + logPZ(Z,alpha);
    %samples = repmat(samples,nIter,1);
    
    %############## RUN GIBBS SAMPLER ###############
    
    for i = 2:nIter
        disp(num2str(i))
        
        %--------- sample alpha -----------------------%
        K = sum(sum(Z)>0);
        alpha = gamrnd(1+K,1/(1+Hn));
        
        %----------- sample sx -----------------------%
        sx2 = gamrnd(1,1);
        score2 = lik(X,Z,sx2,sw);
        score1 = lik(X,Z,sx,sw);
        if rand < exp(score2-score1); sx=sx2; score1=score2; end
        
        %----------- sample sw -----------------------%
        sw2 = exprnd(10);
        score2 = lik(X,Z,sx,sw2);
        if rand < exp(score2-score1); sw=sw2; end
        
        for n = randperm(N)
            
            %------- sample missing data -------------%
            if any(na(n,:))
                [L W] = lik(X,Z,sx,sw);
                z = Z(n,sum(Z)>0);
                X(n,na(n,:)) = z*W(:,na(n,:)) + sx*randn(1,sum(na(n,:)));
            end
            
            m = sum(Z(1:N~=n,:));
            K_active = find(m>0);
            
            %-------- sample active features -------------%
            for k = K_active
                Z2 = Z;
                Z2(n,k) = 1; p2 = log(m(k)/N) + lik(X,Z2,sx,sw);
                Z2(n,k) = 0; p = log(1-m(k)/N) + lik(X,Z2,sx,sw);
                if rand < 1/(1+exp(p-p2))
                    Z(n,k) = 1;
                else
                    Z(n,k) = 0;
                end
            end
            
            %----------- sample new features -------------%
            kn = poissrnd(alpha/N);
            K_inactive = find(m==0);
            kn = min(length(K_inactive),kn);
            Z2 = Z; Z2(n,K_inactive(1:kn)) = 1;
            score2 = lik(X,Z2,sx,sw);
            score1 = lik(X,Z,sx,sw);
            if rand < exp(score2-score1)
                Z = Z2; score1 = score2;
            end 
        end
        
        %store feature matrix
        if missing; samples(i).X = X; end
        samples(i).Z = Z;
        samples(i).score = score1 + logPZ(Z,alpha);
        
        if mod(i,10)==0 && ~isempty(savefile)
            save(savefile,'samples');
        end
    end