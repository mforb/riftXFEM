ld = dir('PRES_xmas_tip*');
results_path = './PRES_xmas_PP';
mkdir(results_path);
global results_path
global zoom_dim
global Hidden
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

t = tiledlayout(2,2,'TileSpacing','Compact');
% tile 1
nexttile
plot([1:length(knm1)],knm1)
xlabel('step')
title('SIFs end 1')
legend({'K1','K2'})

nexttile
plot([1:length(knm2)],knm2)
xlabel('step')
title('SIFs end 2')
legend({'K1','K2'})

nexttile
plot([1:length(t1)],t1)
title('propagation angle, end 1')
xlabel('step')

nexttile
plot([1:length(t2)],t2)
title('propagation angle, end 2')
xlabel('step')
%plotMesh(node+dfa*[uxAna uyAna],element,elemType,'r-',plotNode)

figure_name = ['Knum_results'];
print([results_path,'/',figure_name],'-dpng','-r300')

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
  plotFieldXfemT3(xCrk,pos,enrichNode,crackNode,u,...
    elemCrk,vertexElem,cornerElem,splitElem,tipElem,xVertex,xTip,typeElem,1) ;
end

%plots of the first time-step
if 1
  dname = ld(end).name;
  lname = [dname,'/crack5.mat']; 
  load(lname)
  zoom_dim(1,:) = [min(xCrk.coor(:,1))-20000,max(xCrk.coor(:,1))+20000];
  zoom_dim(2,:) = [min(xCrk.coor(:,2))-10000,max(xCrk.coor(:,2))+10000];
  plotFieldXfemT3(xCrk,pos,enrichNode,crackNode,u,...
    elemCrk,vertexElem,cornerElem,splitElem,tipElem,xVertex,xTip,typeElem,66) ;
end

  


