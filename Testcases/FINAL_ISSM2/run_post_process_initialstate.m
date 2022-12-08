% which run do you want to process 
ld = dir('./FINAL/PRES_xmas_tip1_10*');
% where do you want it stored
results_path = './FINAL/PRES_PP';
mkdir(results_path);

path(path,'../../')
path(path,'../../Crackprocessing')
path(path,'../../Mesh')
path(path,'../../Routines_XFEM')
path(path,'../../Routines_ICEM')
path(path,genpath('~/Softs/MATLAB/TOOLS/'));
fontSize1 = 14; 
fontSize2 = 12; 
mag       = 4000;
plotfields = 1;

global results_path
global zoom_dim
global Hidden
global fontSize2 fontSize1
global elemType 
Hidden = 0;
global E C Cm1 nu P
global melange melangeforce wall_int epsilon
epsilon = 5;
wall_int = 2
epsilon = 5
melange = 0
melangeforce = 0
E = 9.6e9; nu = 0.33; P = 1 ;
elemType = 'T3';
sigmato = P ;
stressState = 'PlaneStrain' ;
if( strcmp(stressState,'PlaneStress') )
    C = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
else
    C = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
    Cm1 = E*0.1/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
end

% the original crack geometry
srift2 = shaperead('~/Work/Shapefiles/rift_2005.shp');
xs = srift2.X
ys = srift2.Y
xs(end) = []; %get rid of trailin NaN
ys(end) = [];
xCr_original.coor = [fliplr(xs)',fliplr(ys)'] 

tip1 = [ ones(1,8), zeros(1,6), ones(1,2)];
tip2 = [ zeros(1,8), ones(1,6), ones(1,2)];

knm1 = [];
knm2 = [];
t1 =[];
t2 =[];
% read all of the SIF values
for i = 1:length(ld)
  dname = ld(i).name;
  lname = [dname,'/crack.mat']; 
  load(lname)
  knm1 = [knm1,Knumerical{1}];
  knm2 = [knm2,Knumerical{2}];
  t1 = [t1,ThetaInc{1}];
  t2 = [t2,ThetaInc{2}];
end

  zoom_dim(1,:) = [min(xCr.coor(:,1))-3000,max(xCr.coor(:,1))+3000];
  zoom_dim(2,:) = [min(xCr.coor(:,2))-3000,max(xCr.coor(:,2))+3000];


%plots of the first time-step
%reset C
if( strcmp(stressState,'PlaneStress') )
    C = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
else
    C = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
    Cm1 = E*0.1/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
end
%plots of the first time-step
if 1
  dname = [ld(1).name];
  lname = ['./FINAL/',dname,'/crack1.mat']; 
  load(lname)
  TR = triangulation(element,node);
  if plotfields
  [ca,cax,cay] = plotFieldXfemT3_pp(xCrk,pos,enrichNode,crackNode,u,...
    elemCrk,vertexElem,cornerElem,splitElem,tipElem,xVertex,xTip,typeElem,66) ;
  end
  fu = full(u);
  numnode = length(node);
  Stdux = fu(1:2:2*numnode) ;
  Stduy = fu(2:2:2*numnode) ;
  %[crackLips,flagP] = f_cracklips( u, xCrk, enrDomain, typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);

  f = figure();
  f.Position = [ 0, 0, 1200, 700];
  hold on
  % the crack is saved after a propagation step, so we need to modify the crack to plot 
  %xCrk(1).coor(1,:)=[];
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
  if ~isempty(zoom_dim)
    xlim(zoom_dim(1,:));
    ylim(zoom_dim(2,:));
    yticks(-1170000:10000:-1100000);
    xticks(-20000:10000:100000);
    f_publish_fig(f,'s');
    figure_name = ['crackwalls',num2str(mag),'_start_zoom'];
    print([results_path,'/',figure_name],'-dpng')
  end
  clf();
  %f_plot_wall_forces(u,xCrk,[],typeElem,elemForce,elemGap,elemCrk,splitElem,vertexElem,tipElem,1)
  [~,ylg,yls]=f_plot_wf(u,xCrk,[],typeElem,elemForce,elemGap,elemCrk,splitElem,vertexElem,xVertex,tipElem,xTip,crackNode,typeElem,enrichNode,pos,1);
  clf();
  f = figure();
  f.Position = [0 0 1200 700 ];
  trisurf(element,node(:,1),node(:,2),Stduy)
  trisurf(element,node(:,1),node(:,2),Stduy)
  axis equal; view(2); shading interp; cb = colorbar();
  cb.Label.String = "displacement";
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



