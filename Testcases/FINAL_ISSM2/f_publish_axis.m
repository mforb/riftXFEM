function [ b ] = f_publish_axis( ax,ss )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Tue Sep 7 2022
switch(ss)
case 't'
  xlim manual
  ylim([-1.36,-0.99]*1e6);
  xlim([-3.22,3.767]*1e5);
case 's'
end
ax.Box='off';
ax.Title = text();
ax.FontSize = 16;
ax.LineWidth = 1.2;
ax.Color = 'none';
ax.TickDir = 'out';
ax.TickLength = [ 0.005 0.01 ];
% change ticklabels to km
xt = ax.XAxis.TickValues;
xt = xt/1000;
ax.XAxis.TickLabels = strsplit(num2str(xt));
ax.XAxis.Exponent = 0;
% 
yt = ax.YAxis.TickValues;
yt = yt/1000;
ax.YAxis.TickLabels = strsplit(num2str(yt));
ax.YAxis.Exponent = 0;

tpos = plotboxpos(ax);
b = annotation("rectangle",tpos,LineWidth=1.2)
grid(ax,'on');

% 


%b = axes('position',ax.Position,'innerposition',ax.InnerPosition,'box','on','ytick',[],'xtick',[],'color','none','linewidth',1)



