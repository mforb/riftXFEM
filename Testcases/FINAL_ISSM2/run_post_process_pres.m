path(path,'../../')
path(path,'../../Crackprocessing')
path(path,'../../Mesh')
path(path,'../../Routines_XFEM')
path(path,'../../Routines_ICEM')
path(path,genpath('~/Softs/MATLAB/TOOLS/'));
fontSize1 = 14; 
fontSize2 = 12; 
mag       = 2000;

ld = dir('PRES_xmas_tip1_10*');
results_path = './PRES_PP';
mkdir(results_path);
global results_path
global zoom_dim
global Hidden
global fontSize2 fontSize1
global elemType 
Hidden = 0;
global E C nu P
global melange
melange = 0
melangeforce = 0
E = 9.6e9; nu = 0.3; P = 1 ;
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

  zoom_dim(1,:) = [min(xCrk.coor(:,1))-5000,max(xCrk.coor(:,1))+5000];
  zoom_dim(2,:) = [min(xCrk.coor(:,2))-5000,max(xCrk.coor(:,2))+5000];


%plots of the first time-step
if 0
  dname = ld(end).name;
  lname = [dname,'/crack2.mat']; 
  load(lname)
  TR = triangulation(element,node);
  zoom_dim(1,:) = [min(xCrk.coor(:,1))-5000,max(xCrk.coor(:,1))+5000];
  zoom_dim(2,:) = [min(xCrk.coor(:,2))-5000,max(xCrk.coor(:,2))+5000];
  [ca,cax,cay] = plotFieldXfemT3_pp(xCrk,pos,enrichNode,crackNode,u,...
    elemCrk,vertexElem,cornerElem,splitElem,tipElem,xVertex,xTip,typeElem,66) ;
  fu = full(u);
  numnode = length(node);
  Stdux = fu(1:2:2*numnode) ;
  Stduy = fu(2:2:2*numnode) ;

  f = figure();
  hold on
  [crackLips,flagP,elemGap] = f_find_cracklips( u, xCrk, 1, [], typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
  dfac = 1 ;
  triplot(TR);
  hold on
  axis equal;
  f_plotCrack_pp(crackLips,mag)
  print([results_path,'/crackwalls',num2str(mag),'_end'],'-dpng','-r300')
  if ~isempty(zoom_dim)
    xlim(zoom_dim(1,:));
    ylim(zoom_dim(2,:));
    figure_name = ['crackwalls',num2str(mag),'_end_zoom'];
    print([results_path,'/',figure_name],'-dpng','-r300')
  end
  clf();
  f_plot_wall_forces(u,xCrk,[],typeElem,elemForce,elemGap,elemCrk,splitElem,vertexElem,tipElem,66)
  clf();
  trisurf(element,node(:,1),node(:,2),Stduy)
  axis equal; view(2); shading interp; colorbar
  cm = flipud(cbrewer2('RdBu', 256));
  colormap(cm);
  caxis([-4,1]);
  title('Y displacement')
  print([results_path,'/end_ydisp'],'-dpng','-r300')
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
  dname = ld(1).name;
  lname = [dname,'/crack1.mat']; 
  load(lname)
  TR = triangulation(element,node);
  %plotFieldXfemT3_pp(xCrk,pos,enrichNode,crackNode,u,...
  %  elemCrk,vertexElem,cornerElem,splitElem,tipElem,xVertex,xTip,typeElem,1,ca,cax,cay);
  fu = full(u);
  numnode = length(node);
  Stdux = fu(1:2:2*numnode) ;
  Stduy = fu(2:2:2*numnode) ;
  %[crackLips,flagP] = f_cracklips( u, xCrk, enrDomain, typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);

  % the crack is saved after a propagation step, so we need to modify the crack to plot 
  xCrk(1).coor(1,:)=[];

  f = figure();
  hold on
  [crackLips,flagP,elemGap] = f_find_cracklips( u, xCrk, 1, [], typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
  dfac = 1 ;
  triplot(TR);
  hold on
  axis equal;
  f_plotCrack_pp(crackLips,mag)
  print([results_path,'/crackwalls',num2str(mag),'_start'],'-dpng','-r300')
  if ~isempty(zoom_dim)
    xlim(zoom_dim(1,:));
    ylim(zoom_dim(2,:));
    figure_name = ['crackwalls',num2str(mag),'_start_zoom'];
    print([results_path,'/',figure_name],'-dpng','-r300')
  end
  clf();
  f_plot_wall_forces(u,xCrk,[],typeElem,elemForce,elemGap,elemCrk,splitElem,vertexElem,tipElem,1)
  clf();
  trisurf(element,node(:,1),node(:,2),Stduy)
  axis equal; view(2); shading interp; colorbar
  cm = flipud(cbrewer2('RdBu', 256));
  colormap(cm);
  caxis([-4,1]);
  title('Y displacement')
  print([results_path,'/start_ydisp'],'-dpng','-r300')
  clf();
end

  


