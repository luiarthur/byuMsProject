function plot_ibp_matrices(X,Z,Y)
clf;
subplot(1,3,1)
imagesc(X)
v = axis;
title(['X '])
subplot(1,3,2)
imagesc(Z)
% axis(v)
title('Z ')
subplot(1,3,3)
imagesc(Y)
axis(v)
title('Y ')