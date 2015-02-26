
function plot_and_save_nips_graphs(X,Y,Z,N,Zsamples,lP, EZZt, save_figs, save_dir)
if(nargin<7)
    SAVEFIGS =0;
else
    SAVEFIGS = save_figs;
end
figure(1)
imagesc(X)
% title('Training data X')
colormap([0 0 0; 1 1 1]);
set(gca,'YTick',1:N);
set(gca,'FontSize',14);
set(gca,'LineWidth',2');
xlabel('T')
ylabel('N')
axis fill
if(SAVEFIGS)
print(gcf, '-deps2',  [save_dir 'X.eps']);
end
figure(2)
plot_graph(Z)
axis fill
if(SAVEFIGS)
    print(gcf, '-deps2',  [save_dir 'graph.eps']);
end
% title('Disconnected generative graph')
figure(3)
imagesc(Z*Z')
set(gca,'XTick',[])
set(gca,'YTick',[])

colormap('hot')
set(gca,'FontSize',14);
c = colorbar
set(c,'FontSize',14);
axis fill
% title('Z*Z'' training data')
if(SAVEFIGS)

    print(gcf, '-depsc2',  [save_dir 'ZZt-actual.eps']);
end


figure(5)
imagesc(EZZt)
set(gca,'XTick',[])
set(gca,'YTick',[])

colormap('hot')
c = colorbar
set(c,'FontSize',14);
axis fill
% title('Z*Z'' training data')
if(SAVEFIGS)

    print(gcf, '-depsc2',  [save_dir 'EZZt-samples.eps']);
end

figure(7)
plot_graph(Zsamples{min(find(lP==max(lP)))})
axis fill
if(SAVEFIGS)

    print(gcf, '-deps2',  [save_dir 'mlZ.eps']);
end

