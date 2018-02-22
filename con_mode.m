function [cs,h] = con_mode(x,y,mode,maxm,dbsep,dbmax,t);


x = real(x);
y = real(y);
mode = abs(mode');
v = (0:-dbsep:-dbmax)';
[cs,h] = contour(x,y,20*log10(mode/maxm),v);
xlim([min(x),max(x)]);
ylim([min(y),max(y)]);
xlabel('x'); 
ylabel('y'); 
title(t);
% v = axis;
% set(gca,'PlotBoxAspectRatio',[v(2)-v(1) v(4)-v(3) 1]);
