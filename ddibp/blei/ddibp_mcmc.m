function samples = ddibp_mcmc(X,D,f,nIter,sx,sw,alpha,savefile)
    
    % MCMC algorithm for approximate inference in the distance-dependent
    % Indian buffet process (dd-IBP).
    %
    % USAGE: samples = ddibp_mcmc(X,D,f,nIter,sx,sw,alpha,[savefile])
    %
    % INPUTS:
    %   X - [observations x dimensions] data matrix (set missing data to NaN)
    %   D - [observations x observations] distane matrix
    %   f - decay function handle
    %   nIter - number of MCMC iterations
    %   sx - observation noise
    %   sw - prior noise
    %   alpha - concentration parameter
    %   savefile (optional) - file to save intermediate results to
    %
    % OUTPUTS:
    %   samples - [nIter x 1] structure containing the following fields:
    %               .score - complete log-likelihood
    %               .Z - latent feature matrix
    %               .X - complete data matrix (if missing data)
    %
    % Sam Gershman, last updated July 2013
    
    %############ INITIALIZATION #############
    
    N = size(X,1);
    na = isnan(X);
    missing = any(na(:));
    if missing
        X(na) = randn*sx;
    end
    
    %--------- parameters --------------%
    K = size(X,2);
    %K = 100;      %maximal number of features (set this to be large)
    
    %--------- compute normalized proximity matrix -----------%
    S = f(D);   %proximity matrix
    h = zeros(1,N);
    for n = 1:N
        S(n,n) = 1;
        h(n) = sum(S(n,:));
        S(n,:) = S(n,:)./h(n);
    end
    logS = log(S);
    Hn = sum(1./h);
    
    %------- initialize customer assignments ------------%
    disp(num2str(1))
    C = zeros(N,K); Z = zeros(N,K);
    lambda = zeros(N,1);
    Kn = zeros(1,K);    %feature ownership vector
    for n = 1:N
        lambda(n) = poissrnd(alpha/h(n));
        if n==1; lambda(n) = max(lambda(n),1); end % make sure there's at least one active feature
        if lambda(n) > 0
            ix = find(Kn==0);
            ix = ix(1:min(lambda(n),length(ix)));
            ix(ix>K)=[];
            Kn(ix) = n;
        end
        for k = 1:K
            C(n,k) = fastrandsample(S(n,:));
        end
    end
    for k = find(Kn); Z(:,k) = C_to_Z(C(:,k),Kn(k)); end
    
    %pre-allocate
    lp = zeros(1,N);
    
    %---------- initialize results structure -----------%
    samples.accept = 0;
    samples.Z = Z;
    samples.score = lik(X,Z,sx,sw) + logPZ_ddibp(Kn,logS,C,alpha,h);
    %samples = repmat(samples,nIter,1);
    
    %############## RUN GIBBS SAMPLER ###############
    
    for i = 2:nIter
        disp(num2str(i));
        
        %----------- sample alpha --------------------%
        alpha = gamrnd(1+sum(lambda),1/(1+Hn));
        
        %----------- sample sx -----------------------%
        sx2 = gamrnd(1,1);
        score2 = lik(X,Z,sx2,sw);
        score1 = lik(X,Z,sx,sw);
        if rand < exp(score2-score1); sx=sx2; score1 = score2; end
        
        %----------- sample sw -----------------------%
        sw2 = exprnd(10);
        score2 = lik(X,Z,sx,sw2);
        if rand < exp(score2-score1); sw=sw2; end
        
        samples(i).accept = 0;
        for n = randperm(N)
            
            %------- sample missing data -------------%
            if any(na(n,:))
                [L W] = lik(X,Z,sx,sw);
                z = Z(n,sum(Z)>0);
                X(n,na(n,:)) = z*W(:,na(n,:)) + sx*randn(1,sum(na(n,:)));
            end
            
            v = 1:N; v(n) = []; %all objects except n
            
            %-------- sample active features -------------%
            for k = find(Kn==n)     %loop over features owned by object n
                for m = v           %loop over all objects except n
                    
                    %pre-compute likelihood
                    C2 = C; Z2 = Z;
                    C2(m,k)=n; Z2(:,k) = C_to_Z(C2(:,k),n); lik1 = lik(X,Z2,sx,sw);
                    C2(m,k)=m; Z2(:,k) = C_to_Z(C2(:,k),n); lik0 = lik(X,Z2,sx,sw);
                    
                    for j = 1:N    %loop over possible customer assignments
                        if reachable(C(:,k),j,n)
                            lp(j) = logS(m,j) + lik1;
                        else
                            lp(j) = logS(m,j) + lik0;
                        end
                    end
                    p = exp(lp-logsumexp(lp,2));
                    C(m,k) = fastrandsample(p);
                    Z(:,k) = C_to_Z(C(:,k),n);
                end
            end
            
            %----------- sample new features -------------%
            lambda2 = poissrnd(alpha/h(n));
            df = lambda2 - lambda(n);
            Kn2 = Kn; C2 = C; Z2 = Z;
            if df > 0   %activate new features
                ix = find(Kn==0);                %indices of inactive features
                ix = ix(1:min(df,length(ix)));
                Kn2(ix) = n;
                for k = ix;
                    C2(n,k) = fastrandsample(S(n,:));
                    Z2(:,k) = C_to_Z(C2(:,k),n);
                end
            elseif df < 0        %deactive features
                ix = find(Kn==n);
                df2 = min(length(ix),lambda(n)-lambda2);
                ix = ix(end-df2+1:end);
                ix = ix(randperm(length(ix)));   %scramble order
                Kn2(ix) = 0;
                Z2(:,ix) = 0;
            end
            score2 = lik(X,Z2,sx,sw);
            if df==0; score1=score2; else score1 = lik(X,Z,sx,sw); end
            if rand < exp(score2 - score1)
                Z = Z2; C = C2; Kn = Kn2; lambda(n) = lambda2;
                samples(i).accept = samples(i).accept + 1/N;
            end
        end
        
        %store feature matrix
        samples(i).Z = Z;
        if missing; samples(i).X = X; end
        samples(i).score = lik(X,Z,sx,sw); + logPZ_ddibp(Kn,logS,C,alpha,h);
        
        if nargin >= 8 && mod(i,10)==0
            save(savefile,'samples');
        end
        
    end