function [retZ,retY] = cannonize(Z,Y)

i = (size(Z,1)-1):-1:0;
p = 2.^(i);
sv= p*Z;
    
[v,i] = sort(sv,'descend');

retZ = Z(:,i);
retY = Y(i,:);
