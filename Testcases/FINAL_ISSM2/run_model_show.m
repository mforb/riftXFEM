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


figure(1);
clf
t = tiledlayout(1,2,'TileSpacing','Compact');

% tile 1
nexttile
hold on
grid on
triplot(TR);
axis equal;
xlabel('easting (m)','FontSize',fontSize2)
ylabel('northing (m)','FontSize',fontSize2);
title('ISSM mesh','FontSize',fontSize1)

nexttile
hold on
grid on
patch('faces',element,'vertices',node,'facevertexcdata',vonmises);
axis equal;
shading flat;
cm = cbrewer2('BuPu', 256);
colormap(cm);
set(gca,'ColorScale','log');
colorbar;
caxis([1e4,1e6]);
axis equal;
xlabel('easting (m)','FontSize',fontSize2)
ylabel('northing (m)','FontSize',fontSize2);
title('ISSM Vonmises stress','FontSize',fontSize1)


figure_name = ['issm_stress'];
print([results_path,'/',figure_name],'-dpng','-r200')
saveas(t,[results_path,'/',figure_name],'epsc')

figure(2);
clf
t = tiledlayout(2,2,'TileSpacing','Compact');

% tile 1
nexttile
hold on
grid on
triplot(TR);
axis equal;
xlabel('easting (m)','FontSize',fontSize2)
ylabel('northing (m)','FontSize',fontSize2);
title('ISSM mesh','FontSize',fontSize2)

indx = -3.18e5;
indy = -1.02e6;
yl = get(gca,'ylim');
xl = get(gca,'xlim');
b1 = plot([indx,indx],[yl(1),indy]);
b2 = plot([indx,xl(2)],[indy,indy]);
t = text([indx],[indy],'XFEM subdomain')
set(t,'position',tp - [-8000 30000 0])
set(t,'fontweight','bold')
set(t,'fontsize',fontSize2)
bs = [b1,b2];
set(bs,'linewidth',2);
set(bs,'color',[0.2 0.2 0.2 0.8]);
set(bs,'linestyle','--');
xlim(xl);
ylim(yl);

nexttile
indx = find(cpos(:,1)>-3.18e5);
indy = find(cpos(:,2)<-1.02e6);
in = intersect(indx,indy);
element = element(in,:);
cpos    = cpos(in,:);
TR = triangulation(element,node);
cpos = TR.incenter;
%plotMesh(node,element,elemType,'b-','yes',f)
triplot(TR)
axis equal;
hold on
xlabel('easting (m)','FontSize',fontSize2)
ylabel('northing (m)','FontSize',fontSize2);
title('XFEM starting mesh','FontSize',fontSize2)

cp1 = plot(xCr.coor(:,1),xCr.coor(:,2),'r');
all_node_num = unique(element);
numnode = size(all_node_num,1); 
node = node(all_node_num,:);
for i=1:3*length(element)
  element(i) = find(all_node_num==element(i));
end
edges_front2 = [];
for i=1:length(edges_front)
  if find(all_node_num==element(i,1)) & find(all_node_num==element(i,2));
    edges_front2 = [edges_front2 ; find(all_node_num==edges_front(i)),find(all_node_num==edges_front(i,2))];
  end
end
edges_front = edges_front2;
bc_fix2 = [];
for i = 1:length(bc_fix)
  if find(all_node_num==bc_fix(i))
  bc_fix2 = [ bc_fix2; find(all_node_num==bc_fix(i))];
  end
end
bc_fix = bc_fix2;

TR = triangulation(element,node);
element1 = element;
node1 = node;
cpos = TR.incenter;
% refinement using ameshreF
in = f_find_points_xCr(cpos,xCr,80000);
%keyboard
in1 = in;
%indx = find(cpos(:,1)>-20e3 & cpos(:,1)<110e3 );
%indy = find(cpos(:,2)>-1180e3 & cpos(:,2)<-1080e3 );

path(path,'~/Softs/ameshref/refinement/')
clear TrefineRG 
[node,element] = TrefineRG(node,element,in);

TR = triangulation(element,node);
%keyboard
cpos = TR.incenter;
in = f_find_points_xCr(cpos,xCr,25000);


%indx = find(cpos(:,1)>0e3 & cpos(:,1)<90e3 );
%indy = find(cpos(:,2)>-1160e3 & cpos(:,2)<-1100e3 );

%in = intersect(indx,indy);
[node,element] = TrefineRG(node,element,in);

TR = triangulation(element,node);
cpos = TR.incenter;
in = f_find_points_xCr(cpos,xCr,20000);


%indx = find(cpos(:,1)>10e3 & cpos(:,1)<80e3 );
%indy = find(cpos(:,2)>-1150e3 & cpos(:,2)<-1110e3 );

%in = intersect(indx,indy);
[node,element] = TrefineRG(node,element,in);

TR = triangulation(element,node);

cpos = TR.incenter;
in = f_find_points_xCr(cpos,xCr,5000,14000);


%indx = find(cpos(:,1)>10e3 & cpos(:,1)<80e3 );
%indy = find(cpos(:,2)>-1150e3 & cpos(:,2)<-1110e3 );

%in = intersect(indx,indy);
[node,element] = TrefineRG(node,element,in);

TR = triangulation(element,node);

cpos = TR.incenter;
in = f_find_points_xCr(cpos,xCr,3000,12000);


%indx = find(cpos(:,1)>10e3 & cpos(:,1)<80e3 );
%indy = find(cpos(:,2)>-1150e3 & cpos(:,2)<-1110e3 );

%in = intersect(indx,indy);
[node,element] = TrefineRG(node,element,in);

nexttile
TR = triangulation(element,node);
triplot(TR);
hold on
cp2 = plot(xCr.coor(:,1),xCr.coor(:,2),'r');
axis equal
xlabel('easting (m)','FontSize',fontSize2)
ylabel('northing (m)','FontSize',fontSize2);
title('XFEM refined','FontSize',fontSize2)

nexttile
TR = triangulation(element,node);
triplot(TR);
hold on
cp3 = plot(xCr.coor(:,1),xCr.coor(:,2),'r');
axis equal
xlim([min(xCr.coor(:,1))-20000,max(xCr.coor(:,1))+20000])
ylim([min(xCr.coor(:,2))-20000,max(xCr.coor(:,2))+20000])
xlabel('easting (m)','FontSize',fontSize2)
ylabel('northing (m)','FontSize',fontSize2);
title('XFEM refined zoom','FontSize',fontSize2)

cpos = TR.incenter;

%in = f_find_points_xCr(cpos,xCr,4000,6000)
%%indx = find(cpos(:,1)>39e3 & cpos(:,1)<65e3 );
%%indy = find(cpos(:,2)>-1137e3 & cpos(:,2)<-1116e3 );

%%in = intersect(indx,indy);
%[node,element] = TrefineRG(node,element,in);
cp = [cp1,cp2,cp3];
set(cp,'linewidth',2);
set(cp,'color',[0.9 0.1 0.1 0.8]);

print([results_path,'/mesh_xfem'],'-dpng','-r200')




