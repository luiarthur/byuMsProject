alpha = 2.2;
lambda = .1;
epsilon = .1;
p = .1;
T = 50;
tempK = 100;
N = 20;

rand('state',0);


Z = zeros(N,tempK);
X = zeros(N,T);


nfirst = poissrnd(alpha,1);
Z(1,1:nfirst) = 1;
m = sum(Z);
max_c = nfirst;

for(i=2:N)
    rns = rand(1,tempK);
      pr = m/i;
      v = pr>rns;
      Z(i,:) = v;
%       Z(i,1:max_c) = v(1:max_c);

      k_more = poissrnd(alpha/i,1);
      if(k_more >=1)
        Z(i,max_c:(max_c+k_more)) = 1;
      end
      max_c = max_c+k_more;
  m = sum(Z);
end

K = max_c;
Z = Z(:,1:K);

Y = rand(K,T);
zi = find(Y<(1-p));
Y = ones(K,T);
Y(zi)=0;


pX = 1-(1-lambda).^(Z*Y)*(1-epsilon);
flips = rand(size(pX));
X=zeros(size(pX));
X(find(flips<pX))=1;

[Z,Y] = cannonize(Z,Y);
figure(2)
imagesc(X)
title(['X Test Data log(P(X,Y,Z)) = ' num2str(logPXYZ(X,Y,Z,alpha,epsilon,lambda,p))])
figure(3)
imagesc(Z)
title('Z Test Data')
figure(4)
imagesc(Y)
title('Y Test Data')
