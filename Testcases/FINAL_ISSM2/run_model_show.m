path(path,'../../')
path(path,'../../Crackprocessing')
path(path,'../../Mesh')
path(path,'../../Routines_XFEM')
path(path,'../../Routines_ICEM')
path(path,genpath('~/Softs/MATLAB/TOOLS/'));
fontSize1 = 14; 
fontSize2 = 12; 
mag       = 2000;

results_path = './MESH_PP';
mkdir(results_path);
global results_path
global zoom_dim
global Hidden
global fontSize2 fontSize1
global elemType 
Hidden = 0;
global E C nu P
global melange

srift2 = shaperead('~/Work/Shapefiles/rift_2005.shp');
xs = srift2.X
ys = srift2.Y
xs(end) = []; %get rid of trailin NaN
ys(end) = [];
xCr.coor = [fliplr(xs)',fliplr(ys)'] 

load ../FINAL_ISSM/import_issm_holly1
vv = sqrt(Vx.*Vx + Vy.*Vy);
vve = mean(vv(element),2);
%element = element(1:5,:);
% we are going to use triangulation to create a Fintx...
% this is necessary if we plan on doing some refinement in the vicinity of the rift(s)
TR = triangulation(element,node);

cpos = TR.incenter;
stress = [ISSM_xx,ISSM_yy,ISSM_xy];
vonmises  = sqrt( ISSM_xx.^2 +ISSM_yy.^2 -ISSM_xx.*ISSM_yy + 3*ISSM_xy.^2 );

rg(1) = min(vonmises);
rg(2) = max(vonmises);
ormin = floor( log10(abs(rg(1))));
ormax = floor( log10(abs(rg(2))));
mr = max([ormin,ormax]) - 1;
lb = floor(rg(1)/(10^mr)) * 10 ^mr;
ub = ceil(rg(2)/(10^mr)) * 10 ^mr;


f = figure();
f.Position = [0, 0, 900, 800 ]
clf

% tile 1
hold on
grid on
triplot(TR);
axis equal;
yl = get(gca,'ylim')
ylim([yl(1)-10000,yl(2)])
xlabel('Easting (km)','FontSize',fontSize2)
ylabel('Northing (km)','FontSize',fontSize2);
%title('ISSM mesh','FontSize',fontSize1)
indx = -3.18e5;
indy = -1.02e6;
yl = get(gca,'ylim');
xl = get(gca,'xlim');
b1 = plot([indx,indx],[yl(1),indy]);
b2 = plot([indx,xl(2)],[indy,indy]);
% text for cut domain
t = text([indx],[indy],'XFEM subdomain')
tp = t.Position;
set(t,'position',tp - [-10000 310000 0])
%set(t,'fontweight','bold')
set(t,'fontsize',16)
t.Color = [0.2 0.2 0.2 1];

% lines for cut domain
bs = [b1,b2];
set(bs,'linewidth',2);
set(bs,'color',[0.2 0.2 0.2 1]);
set(bs,'linestyle','--');
xlim(xl);
ylim(yl);

% shelf front
[v,so] = sort(node(bc_front,1));
pa = plot(node(bc_front(so(17:end)),1),node(bc_front(so(17:end)),2));
pa.Color = [0.6, 0.15, 0.3];
pa.LineWidth = 1.5;

t2 = text(-250*1e3,-1245*1e3,'ice front')
t2.FontSize = 16;
t2.Color = [0.6, 0.15, 0.3];
t2.Rotation = -14;


tvv = patch('faces',element,'vertices',node,'facevertexcdata',vv,'facecolor','interp');
cm1 = cbrewer2('PuRd', 256);
colormap(cm1);
cb1 = colorbar();
cb1.FontSize = 16;
cb1.Label.String = 'Velocity';
cb1.TickDirection = 'out';

tst = patch('faces',element,'vertices',node,'facevertexcdata',vonmises);
shading flat;
cm2 = cbrewer2('BuPu', 256);
%caxis([1e4,1e6]);

set(tst,'visible','off');
set(tvv,'visible','off');
set(cb1,'visible','off');

xlim(xl);
ylim(yl);

f_publish_fig(f,'b');
figure_name = ['issm_model'];
print([results_path,'/',figure_name],'-dpng','-r300')
saveas(f,[results_path,'/',figure_name],'epsc')
set(t,'visible','off');
%set(t2,'visible','off');
set(bs,'visible','off');
figure_name = ['issm_model2'];
print([results_path,'/',figure_name],'-dpng','-r300')
saveas(f,[results_path,'/',figure_name],'epsc')
set(tvv,'visible','on');
set(cb1,'visible','on');
tvv.FaceColor = 'interp';
figure_name = ['issm_velocity'];
print([results_path,'/',figure_name],'-dpng','-r300')
saveas(f,[results_path,'/',figure_name],'epsc')
set(tvv,'visible','off');
set(cb1,'visible','off');
set(tst,'visible','on');
cb2 = colorbar();
cb2.FontSize = 16;
cb2.Label.String = 'von Mises stress';
cb2.TickDirection = 'out';
set(gca,'ColorScale','log');
colormap(cm2);
figure_name = ['issm_vonmises'];
print([results_path,'/',figure_name],'-dpng','-r300')
saveas(f,[results_path,'/',figure_name],'epsc')
