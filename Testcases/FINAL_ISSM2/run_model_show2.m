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


%f = figure('visible','off');
f = figure();
f.Position = [0 0 1200 700 ]
% tile 1
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
grid on
hold on
xlabel('Easting (km)','FontSize',fontSize2)
ylabel('Northing (km)','FontSize',fontSize2);
%title('XFEM starting mesh','FontSize',fontSize2)

cp1 = plot(xCr.coor(:,1),xCr.coor(:,2),'r');
cp1.LineWidth = 2;
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
yticks(-1300000:100000:-1000000);
f_publish_fig(f,'t');
%print([results_path,'/mesh_xfem1'],'-dpng','-r300')
print([results_path,'/mesh_xfem1'],'-dpng')

clf();
TR = triangulation(element,node);
triplot(TR);
grid on
hold on
cp2 = plot(xCr.coor(:,1),xCr.coor(:,2),'r');
cp2.LineWidth = 2;
axis equal
xlabel('Easting (km)','FontSize',fontSize2)
ylabel('Northing (km)','FontSize',fontSize2);
yticks(-1300000:100000:-1000000);
f_publish_fig(f,'t');
%print([results_path,'/mesh_xfem2'],'-dpng','-r300')
print([results_path,'/mesh_xfem2'],'-dpng')
%title('XFEM refined','FontSize',fontSize2)

clf();
TR = triangulation(element,node);
triplot(TR);
grid on
hold on
cp3 = plot(xCr.coor(:,1),xCr.coor(:,2),'r');
cp3.LineWidth = 2;
axis equal
xlim([min(xCr.coor(:,1))-25000,max(xCr.coor(:,1))+25000])
ylim([min(xCr.coor(:,2))-20000,max(xCr.coor(:,2))+20000])
xlabel('Easting (km)','FontSize',fontSize2)
ylabel('Northing (km)','FontSize',fontSize2);
%{title('XFEM refined zoom','FontSize',fontSize2)%}
yticks(-1150000:10000:-1120000);
f_publish_fig(f,'t');
%print([results_path,'/mesh_xfem3'],'-dpng','-r300')
print([results_path,'/mesh_xfem3'],'-dpng')




