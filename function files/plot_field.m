function plot_field(fnbr, titl_str, x, y, z, plotcircle)
%plots the (x, y, z) values with 
%f = figure number fnbr
%title = titl_str
%%%
%plotcircle is the z-value of the unit circle if you want to
% plot it, and 0 if not

f = figure(fnbr);
clf(f);
title(titl_str)
hold on;
s_interp = surf(x, y, z);
%[M, c] = contour(x, y, z, [0, 0], 'black');
s_interp.EdgeColor = 'none';

%sc = scatter(xpos, ypos, 'black', 'filled');
colorbar;
%caxis([-1, 1]);
%caxis([-10, 10]);


if plotcircle > 0
   circle3d(0,0,plotcircle) 
end
hold off;
end


function h = circle3d(x0,y0,r)
    theta=-pi:0.01:pi;
    x=r*cos(theta);
    y=r*sin(theta);
    plot3(x-x0,y-y0,ones(1,numel(x)), 'black')
end