function plot_graph(Z,xlabels,hightlight_causes)

if(nargin < 2)
    xlabels = cell(0);
end

if(nargin<3)
    hightlight_causes = [];
end

K = size(Z,2);
N = size(Z,1);

of = K/2;
upx = linspace(-of,of,K);

of = N/2;
dpx = linspace(-of,of,N);

cla
hold on
mZ = max(max(Z));
% if(length(xlabels)>0)
%     set(gca,'YLim',[-.5 1])
% end
radius =.25;
for(d = 1:N)
    for(u = [1:K hightlight_causes])

        if(Z(d,u) > 0)
            
            x1 = dpx(d);
            x2 = upx(u);
            y1 = 0+radius;
            y2 = 5-radius;
            r = .05; % the problem is that we don't know what the radius of the glyph is 
            
            theta = atan2((y2-y1),(x2-x1));
            D = sqrt((y2-y1)^2+(x2-x1)^2);
            dn = D -r;
            y3 = dn*sin(theta);
            x3 = dn*cos(theta);
            h = line([x1 x2],[y1 y2]);
%             h = line([upx(u) dpx(d)],[1 0]);
            lss = find(hightlight_causes == u);
            if(isempty(lss))
                lss = -1;
            end
            switch(lss)
                case -1

                    set(h,'Color',[.5 .5 .5]);
                    set(h,'LineWidth',1*(Z(d,u)/mZ));
                case 1
                    set(h,'Color',[0 0 0]);
                    set(h,'LineWidth',2*(Z(d,u)/mZ));
                    set(h,'LineStyle','--');
                case 2
                    set(h,'Color',[.4 .4 .4]);
                    set(h,'LineWidth',2*(Z(d,u)/mZ));
                    set(h,'LineStyle','-');
            end


        end
    end
    if(length(xlabels)>0)
        th = text(dpx(d),-(radius+.1),xlabels(d));
        set(th,'Rotation',90)
        set(th,'HorizontalAlignment','right')
    end
end

for(u = 1:K)
     plot_circle(upx(u),5,radius,1,[0 0 0], [1 1 1]);
end
for(d = 1:N)
     plot_circle(dpx(d),0,radius,1,[0 0 0], [0.5 0.5 0.5]);
end


hold off
set(gca,'YTick',[])
set(gca,'XTick',[])
set(gca,'Box','off')

%
% textx = dpx(end)+(dpx(end)-dpx(end-1));
%
% text(textx,0,'Observed Units')
% text(textx,1,'Hidden Units')