function [ ] = f_publish_axis( ax, f )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Tue Sep 7 2022
ax.Box='off';
ax.Title = text();
ax.FontSize = 16;
ax.LineWidth = 1.2;
ax.Color = 'none';
ax.TickDir = 'out';
ax.TickLength = [ 0.005 0.01 ]
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

b = copyobj(ax,f);
b.XTick = [];
b.YTick = [];
b.Box = 'on';
delete(b.XAxis.Label);
delete(b.YAxis.Label);
delete(b.Children);
%b.InnerPosition = ax.InnerPosition;
b.Position = ax.Position;

% 


%b = axes('position',ax.Position,'innerposition',ax.InnerPosition,'box','on','ytick',[],'xtick',[],'color','none','linewidth',1)



