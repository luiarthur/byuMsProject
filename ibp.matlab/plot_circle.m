function plot_circle(x,y, radius, line_width, line_color, fill_color)
if(nargin < 4)
    line_color = [0 0 0];
end
cy =sin(0:.25:2*pi)*radius+y;
cx = cos(0:.25:2*pi)*radius+x;

if(nargin>5)
fh = fill(cx,cy,fill_color);
end

lh = line(cx,cy);
if(nargin<5)
    set(lh,'LineWidth',2);
else
set(lh,'LineWidth',line_width);
end
set(lh,'Color',line_color);