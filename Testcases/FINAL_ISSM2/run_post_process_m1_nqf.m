path(path,'../../')
path(path,'../../Crackprocessing')
path(path,'../../Mesh')
path(path,'../../Routines_XFEM')
path(path,'../../Routines_ICEM')
path(path,'../../Postprocess')
path(path,genpath('~/Softs/MATLAB/TOOLS/'));
fontSize1 = 14; 
fontSize2 = 12; 
mag       = 4000;
plotfields = 1;

ld = dir('./PLSTRESS/MEL1_NQF_tip*');
pre = ('./PLSTRESS/');
results_path = './PLSTRESS/PPNQF1';
mkdir(results_path);
global results_path
global zoom_dim
global Hidden
global fontSize2 fontSize1
global elemType 
Hidden = 0;
global E C nu P
global melange melangeforce wall_int epsilon
epsilon = 5;
wall_int = 2
melange = 1
melangeforce = 0
E = 9.6e9; nu = 0.33; P = 1 ;
elemType = 'T3';
sigmato = P ;
stressState = 'PlaneStress' ;
if( strcmp(stressState,'PlaneStress') )
    C = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
    Cm1 =0.1*E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
else
    C = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
    Cm1 = E*0.1/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
end

% the original crack geometry
srift2 = shaperead('./o_rift.shp');
xs = srift2.X
ys = srift2.Y
xs(end) = []; %get rid of trailin NaN
ys(end) = [];
xCr_original.coor = [xs',ys'] 
xCr_original.melange = ones(length(xCr_original.coor)-1,1);
%xCr(1).melange(1) = 0;
%xCr(1).melange(end) = 0;
xCr_original.width = [0 10 60 150 200 80 0 ] ;

tip1 = [ ones(1,16), zeros(1,12), ones(1,4)];
tip2 = [ zeros(1,16), ones(1,12), ones(1,4)];

knm1 = [];
knm2 = [];
t1 =[];
t2 =[];
% read all of the SIF values
for i = 1:length(ld)
  dname = ld(i).name;
  lname = [pre,dname,'/crack.mat']; 
  load(lname)
  knm1 = [knm1,Knumerical{1}];
  knm2 = [knm2,Knumerical{2}];
  t1 = [t1,ThetaInc{1}];
  t2 = [t2,ThetaInc{2}];
end

c1 = cbrewer2('set2',4);
c2 = cbrewer2('dark2',4);
rg(1) = min([knm1(1,:),knm1(2,:),knm2(1,:),knm2(2,:)]);
rg(2) = max([knm1(1,:),knm1(2,:),knm2(1,:),knm2(2,:)]);
ormin = floor( log10(abs(rg(1))));
ormax = floor( log10(abs(rg(2))));
mr = max([ormin,ormax]) - 1;
lb = floor(rg(1)/(10^mr)) * 10 ^mr;
ub = ceil(rg(2)/(10^mr)) * 10 ^mr;

t1_cu = cumsum(t1.*tip1);
t2_cu = cumsum(t2.*tip2);

f = figure();
t = tiledlayout(2,2,'TileSpacing','Compact');
f.Position = [ 0, 0, 1200, 900 ];

% tile 1
nexttile
hold on
grid on
for i = 1:2
  plot([1:length(knm1)],knm1(i,:)/1e6,'color',c1(i,:),'linewidth',3)
end
plot([16,16],[lb,ub]/1e6,'color',[30,30,30,200]/255,'linewidth',1)
plot([28,28],[lb,ub]/1e6,'color',[30,30,30,200]/255,'linewidth',1)
ylim([lb,ub]/1e6);
xlim([1,length(knm2)]);
xlabel('Step')
ylabel(['SIF ($\frac{MPa}{\sqrt{m}}$)'],'interpreter','latex','FontSize',14)
%title('SIFs','FontSize',fontSize1)
l1 = legend({'K1','K2'})
ax = gca();
ax.FontSize = 14;

nexttile
hold on
grid on
for i = 1:2
  plot([1:length(knm2)],knm2(i,:)/1e6,'color',c2(i,:),'linewidth',3)
end
plot([16,16],[lb,ub]/1e6,'color',[30,30,30,200]/255,'linewidth',1)
plot([28,28],[lb,ub]/1e6,'color',[30,30,30,200]/255,'linewidth',1)
ylim([lb,ub]/1e6);
xlim([1,length(knm2)]);
xlabel('Step')
ylabel(['SIF ($\frac{MPa}{\sqrt{m}}$)'],'interpreter','latex','FontSize',14)
%title('SIFs','FontSize',fontSize1)
l2 = legend({'K1','K2'})
ax = gca();
ax.FontSize = 14;

nexttile
hold on
plot([1:length(t1_cu)],t1_cu,'color',c1(4,:),'linewidth',3)
plot([1:length(t1)],t1,'color',c1(3,:),'linewidth',3)
plot([16,16],[-pi/3,pi/3],'color',[30,30,30,200]/255,'linewidth',1)
plot([28,28],[-pi/3,pi/3],'color',[30,30,30,200]/255,'linewidth',1)
grid on
ylim([-pi/3,pi/3]);
xlim([1,length(knm2)]);
xlabel('Step');
ylabel('Angle ($^{\circ}$)','interpreter','latex','FontSize',14);
%title('propagation angle','FontSize',fontSize1)
legend({'cumulative angle','angle'})
ax = gca();
ax.FontSize = 14;

nexttile
hold on
plot([1:length(t2_cu)],t2_cu,'color',c2(4,:),'linewidth',3)
plot([1:length(t2)],t2,'color',c2(3,:),'linewidth',3)
plot([16,16],[-pi/3,pi/3],'color',[30,30,30,200]/255,'linewidth',1)
plot([28,28],[-pi/3,pi/3],'color',[30,30,30,200]/255,'linewidth',1)
grid on
ylim([-pi/3,pi/3]);
%title('propagation angle','FontSize',fontSize1)
xlim([1,length(knm2)]);
xlabel('Step');
ylabel('Angle ($^{\circ}$)','interpreter','latex','FontSize',14);
%title('propagation angle','FontSize',fontSize1)
legend({'cumulative angle','angle'})
%plotMesh(node+dfa*[uxAna uyAna],element,elemType,'r-',plotNode)
ax = gca();
ax.FontSize = 14;
%keyboard
%set(l1,'Position',l1.Position - [0, 0.1 0 0 ])
%set(l2,'Position',l2.Position - [0, 0.1 0 0 ])




figure_name = ['Knum_results_M1'];
print([results_path,'/',figure_name],'-dpng','-r300')
saveas(t,[results_path,'/',figure_name],'epsc')


figure_name = ['Knum_results_M1'];
print([results_path,'/',figure_name],'-dpng','-r300')
saveas(t,[results_path,'/',figure_name],'epsc')

%in the last file we loaded the final crack geometry
xCr_final = xCr;
shapefile_name = 'final_rift_m1';
f_crack_shapefile(xCr_final,results_path,shapefile_name);
shapefile_name = 'originalrift_m';
f_crack_shapefile(xCr_original,results_path,shapefile_name);


if 1
  dname = [pre,ld(end).name];
  lname = [dname,'/crack4.mat']; 
  load(lname)
  TR = triangulation(element,node);
  zoom_dim(1,:) = [min(xCrk.coor(:,1))-30000,max(xCrk.coor(:,1))+30000];
  zoom_dim(2,:) = [min(xCrk.coor(:,2))-30000,max(xCrk.coor(:,2))+30000];
  zoom_dim2(1,:) = [min(xCrk.coor(:,1))-3000,max(xCrk.coor(:,1))+3000];
  zoom_dim2(2,:) = [min(xCrk.coor(:,2))-3000,max(xCrk.coor(:,2))+3000];
  if plotfields
  [ca,cax,cay] = plotFieldXfemT3_pp(xCrk,pos,enrichNode,crackNode,u,...
    elemCrk,vertexElem,cornerElem,splitElem,tipElem,xVertex,xTip,typeElem,66) ;
  end
  fu = full(u);
  numnode = length(node);
  Stdux = fu(1:2:2*numnode) ;
  Stduy = fu(2:2:2*numnode) ;

  % the crack is saved after a propagation step, so we need to modify the crack to plot 
  xCrk(1).coor(1,:)=[];
  xCrk(1).coor(end,:)=[];
  f = figure();
  f.Position = [0 0 1200 700 ]
  hold on
  [crackLips,flagP,elemGap] = f_find_cracklips( u, xCrk, 1, [], typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
  dfac = 1 ;
  %triplot(TR);
  hold on
  axis equal;
  f_plotCrack_pp(crackLips,mag,xCr_original)
  ylabel('Northing (km)');
  xlabel('Easting (km)');
  ax = gca();
  ax.FontSize = 16;
  b = f_publish_fig(f,'t');
  print([results_path,'/crackwalls',num2str(mag),'_end'],'-dpng')
  delete(b);
  if ~isempty(zoom_dim2)
    xlim(zoom_dim2(1,:));
    ylim(zoom_dim2(2,:));
    yticks(-1170000:10000:-1100000);
    xticks(-20000:10000:100000);
    f_publish_fig(f,'s');
    figure_name = ['crackwalls',num2str(mag),'_end_zoom'];
    print([results_path,'/',figure_name],'-dpng')
  end
  clf();
  [~,~,ylg,yls]=f_plot_wf(u,xCrk,[],typeElem,elemForce,elemGap,elemCrk,splitElem,vertexElem,xVertex,tipElem,xTip,crackNode,enrichNode,pos,66);
  clf();
  f = figure();
  f.Position = [0 0 1200 500 ]
  trisurf(element,node(:,1),node(:,2),Stduy)
  axis equal; view(2); shading interp; cb = colorbar();
  cb.Label.String = "displacement (m)";
  cb.FontSize = 16;
  ax = gca();
  ax.FontSize = 16;
  
  ylabel('Northing (km)');
  xlabel('Easting (km)');
  cm = flipud(cbrewer2('RdBu', 256));
  colormap(cm);
  caxis([-4,1]);
  %title('Y displacement')
  yticks(-1300000:100000:-1000000);
  f_publish_fig(f,'t');
  print([results_path,'/end_ydisp'],'-dpng')
  clf();
end

%reset C
if( strcmp(stressState,'PlaneStress') )
    C = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
else
    C = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
    Cm1 = E*0.1/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
end
%plots of the first time-step
if 1
  dname = [pre,ld(1).name];
  lname = [dname,'/crack1.mat']; 
  load(lname)
  TR = triangulation(element,node);
  if plotfields
  plotFieldXfemT3_pp(xCrk,pos,enrichNode,crackNode,u,...
    elemCrk,vertexElem,cornerElem,splitElem,tipElem,xVertex,xTip,typeElem,1,ca,cax,cay);
  end
  fu = full(u);
  numnode = length(node);
  Stdux = fu(1:2:2*numnode) ;
  Stduy = fu(2:2:2*numnode) ;

  f = figure();
  f.Position = [ 0, 0, 1200, 700];
  hold on
  % the crack is saved after a propagation step, so we need to modify the crack to plot 
  xCrk(1).coor(1,:)=[];
  [crackLips,flagP,elemGap] = f_find_cracklips( u, xCrk, 1, [], typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
  dfac = 1 ;
  %triplot(TR);
  hold on
  axis equal;
  f_plotCrack_pp(crackLips,mag,xCr_original)
  ylabel('Northing (km)');
  xlabel('Easting (km)');
  ax = gca();
  ax.FontSize = 16;
  b = f_publish_fig(f,'t');
  print([results_path,'/crackwalls',num2str(mag),'_start'],'-dpng','-r300')
  delete(b);
  if ~isempty(zoom_dim2)
    xlim(zoom_dim2(1,:));
    ylim(zoom_dim2(2,:));
    yticks(-1170000:10000:-1100000);
    xticks(-20000:10000:100000);
    f_publish_fig(f,'s');
    figure_name = ['crackwalls',num2str(mag),'_start_zoom'];
    print([results_path,'/',figure_name],'-dpng')
  end
  clf();
  [~,ylg,yls] = f_plot_wf(u,xCrk,[],typeElem,elemForce,elemGap,elemCrk,splitElem,vertexElem,xVertex,tipElem,xTip,crackNode,enrichNode,pos,1,ylg,yls)
  clf();
  f = figure();
  f.Position = [0 0 1200 700 ];
  trisurf(element,node(:,1),node(:,2),Stduy)
  trisurf(element,node(:,1),node(:,2),Stduy)
  axis equal; view(2); shading interp; cb = colorbar();
  cb.Label.String = "displacement (m)";
  cm = flipud(cbrewer2('RdBu', 256));
  colormap(cm);
  cb.FontSize = 16;
  ax = gca();
  ax.FontSize = 16;
  ylabel('Northing (km)');
  xlabel('Easting (km)');
  caxis([-4,1]);
  %title('Y displacement')
  yticks(-1300000:100000:-1000000);
  f_publish_fig(f,'t');
  print([results_path,'/start_ydisp'],'-dpng')
  clf();
end
