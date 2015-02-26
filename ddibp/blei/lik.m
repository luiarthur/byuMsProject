function [score W] = lik(X,Z,sx,sw)
    %evaluates likelihood
    
    [N M] = size(X);
    Z(:,sum(Z)==0) = [];
    K = size(Z,2);
    H = Z'*Z + (sx/sw)*eye(K);
    dH = det(H); if dH <=0; dH = 0.000001; end
    Q = H\Z';
    score = -0.5*(1/sx)*trace(X'*(eye(N)-Z*Q)*X) - N*M*0.5*log(2*pi);
    score = score - M*(N-K)*log(sqrt(sx)) - K*M*log(sqrt(sw)) - M*0.5*log(dH);
    
    if nargout == 2
        W = Q*X;
    end