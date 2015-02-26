function [X,Z,Y] = generate_test_data(N,T,alpha,lambda,epsilon,p,K_specified)
% [X,Z,Y] = generate_test_data(N,T,alpha,lambda,epsilon,p,K_specified)
% generates noisy-or test data with given model parameters
% default values for parameters are K_specified=-1;, p = .2; 
% epsilon = .0001;lambda = .5;alpha = 2.2;T = 1000;N = 20;
% if K_specified is positive then Z and Y will have hidden dimension K
% otherwise K will be random
%
debug = 0;
switch(nargin)
    case 0
        K_specified=-1;
        p = .2;
        epsilon = .0001;
        lambda = .5;
        alpha = 2.2;
        T = 1000;
        N = 20;
    case 1
        K_specified=-1;
        p = .2;
        epsilon = .0001;
        lambda = .5;
        alpha = 2.2;
        T = 1000;

    case 2
        K_specified=-1;
        p = .2;
        epsilon = .0001;
        lambda = .5;
        alpha = 2.2;

    case 3
        K_specified=-1;
        p = .2;
        epsilon = .0001;
        lambda = .5;

    case 4
        K_specified=-1;
        p = .2;
        epsilon = .0001;

    case 5
        K_specified=-1;
        p = .2;
    case 6
        K_specified=-1;

end

[Z,K] = ibp_generate(N,alpha,K_specified);
X = zeros(N,T);

Y = rand(K,T);
zi = find(Y<(1-p));
Y = ones(K,T);
Y(zi)=0;


pX = 1-((1-lambda).^(Z*Y))*(1-epsilon);
flips = rand(size(pX));
X=zeros(size(pX));
X(find(flips<pX))=1;

[Z,Y] = cannonize(Z,Y);
if(debug)
    plot_ibp_matrices(X,Z,Y)
end
