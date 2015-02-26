function [retZ,retY] = clean(Z,Y)

columncounts = sum(Z,1);
nzc = find(columncounts~=0);
retZ = Z(:,nzc);
retY = Y(nzc,:);

