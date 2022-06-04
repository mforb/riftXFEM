path(path,'../../')
path(path,'../../Crackprocessing')
path(path,'../../Mesh')
path(path,'../../Routines_XFEM')
fontSize1 = 14; 
fontSize2 = 12; 

ld = dir('ISSM_xmas_tip*');
results_path = './ISSM_xmas_PP';
mkdir(results_path);
global results_path
global zoom_dim
global Hidden
global fontSize2 fontSize1
Hidden = 0;
global E C nu P
E = 9.6e9; nu = 0.3; P = 1 ;
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


t = tiledlayout(2,2,'TileSpacing','Compact');

% tile 1
nexttile
hold on
grid on
plot([9,9],[lb,ub],'color',[30,30,30,200]/255,'linewidth',1)
plot([17,17],[lb,ub],'color',[30,30,30,200]/255,'linewidth',1)
for i = 1:2
  plot([1:length(knm1)],knm1(i,:),'color',c1(i,:),'linewidth',3)
end
ylim([lb,ub]);
xlim([1,length(knm2)]);
xlabel('step','FontSize',fontSize2)
title('SIFs','FontSize',fontSize1)
l = legend({'K1','K2'})
l.FontSize = fontSize2;

nexttile
hold on
grid on
plot([9,9],[lb,ub],'color',[30,30,30,200]/255,'linewidth',1)
plot([17,17],[lb,ub],'color',[30,30,30,200]/255,'linewidth',1)
for i = 1:2
  plot([1:length(knm2)],knm2(i,:),'color',c2(i,:),'linewidth',3)
end
ylim([lb,ub]);
xlim([1,length(knm2)]);
xlabel('step','FontSize',fontSize2)
title('SIFs','FontSize',fontSize1)
l = legend({'K1','K2'})
l.FontSize = fontSize2;

nexttile
hold on
plot([9,9],[-pi/3,pi/3],'color',[30,30,30,200]/255,'linewidth',1)
plot([17,17],[-pi/3,pi/3],'color',[30,30,30,200]/255,'linewidth',1)
plot([1:length(t1_cu)],t1_cu,'color',c1(4,:),'linewidth',3)
plot([1:length(t1)],t1,'color',c1(3,:),'linewidth',3)
grid on
ylim([-pi/3,pi/3]);
xlim([1,length(knm2)]);
xlabel('step','FontSize',fontSize2)
title('propagation angle','FontSize',fontSize1)
legend({'cumul angle','angle'})
l.FontSize = fontSize2;

nexttile
hold on
plot([9,9],[-pi/3,pi/3],'color',[30,30,30,200]/255,'linewidth',1)
plot([17,17],[-pi/3,pi/3],'color',[30,30,30,200]/255,'linewidth',1)
plot([1:length(t2_cu)],t2_cu,'color',c2(4,:),'linewidth',3)
plot([1:length(t2)],t2,'color',c2(3,:),'linewidth',3)
grid on
ylim([-pi/3,pi/3]);
title('propagation angle','FontSize',fontSize1)
xlim([1,length(knm2)]);
legend({'cumul angle','angle'})
l.FontSize = fontSize2;
%plotMesh(node+dfa*[uxAna uyAna],element,elemType,'r-',plotNode)

figure_name = ['Knum_results'];
print([results_path,'/',figure_name],'-dpng','-r300')
saveas(t,[results_path,'/',figure_name],'epsc')

%in the last file we loaded the final crack geometry
xCr_final = xCr;
srift_final = srift2;
srift_final.BoundingBox = [min(xCr.coor(:,1)), min(xCr.coor(:,2)); max(xCr.coor(:,1)), max(xCr.coor(:,2))];
srift_final.X = xCr_final.coor(:,1)';
srift_final.Y = xCr_final.coor(:,2)';
shapefile_name = 'final_rift';
shapewrite(srift_final,[results_path,'/',shapefile_name]);


%plots of the first time-step
if 1
  dname = ld(1).name;
  lname = [dname,'/crack1.mat']; 
  load(lname)
  zoom_dim(1,:) = [min(xCrk.coor(:,1))-20000,max(xCrk.coor(:,1))+20000];
  zoom_dim(2,:) = [min(xCrk.coor(:,2))-10000,max(xCrk.coor(:,2))+10000];
  [ca,cax,cay] = plotFieldXfemT3_pp(xCrk,pos,enrichNode,crackNode,u,...
    elemCrk,vertexElem,cornerElem,splitElem,tipElem,xVertex,xTip,typeElem,1) ;
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
  dname = ld(end).name;
  lname = [dname,'/crack2.mat']; 
  load(lname)
  zoom_dim(1,:) = [min(xCrk.coor(:,1))-20000,max(xCrk.coor(:,1))+20000];
  zoom_dim(2,:) = [min(xCrk.coor(:,2))-10000,max(xCrk.coor(:,2))+10000];
  plotFieldXfemT3_pp(xCrk,pos,enrichNode,crackNode,u,...
    elemCrk,vertexElem,cornerElem,splitElem,tipElem,xVertex,xTip,typeElem,66,ca) ;
end

  


