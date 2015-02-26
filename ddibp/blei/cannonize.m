function retZ = cannonize(Z)

i = (size(Z,1)-1):-1:0;
p = 2.^(i);
sv= p*Z;
    
[v,i] = sort(sv,'descend');

retZ = Z(:,i);