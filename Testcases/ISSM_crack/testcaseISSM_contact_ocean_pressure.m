%path and variable initialisation
clear all
close all
clc
format long
tic
%---- Define path for subroutines
path(path,'../../')
path(path,'../../Crackprocessing')
path(path,'../../Mesh')
path(path,'../../Routines_XFEM')
path(path,'../../Routines_ICEM')

%define (and make) a path for results
results_path = './ISSM_xmas_soft_stab_nov_superref';
mkdir(results_path);
%copyfile('../Testcase.m',results_path);

%declare global variables here
global L D E nu P C sigmato
global elemType typeMesh typeProblem typeCrack stressState
global xCr xCr_orig deltaInc numstep numcrack
global plotmesh plotNode
global node element numnode numelem bcNodes edgNodes
global results_path fmesh
% global variables for stress
global loadstress FintX FintY FintXY FintH
% global variables for conversion between two coordinate systems
global Rtip QT xTip Tfact
global ISSM_xx ISSM_yy ISSM_xy
global OPT Hidden epsilon melange melangeforce Cm1 xM rift_wall_pressure contact Kpen stabilize
global wall_int stabilize skip_vertex
epsilon = 5 

OPT = 2; Hidden = true;

same_coords = 1

rift_wall_pressure = 1
melange = 0 
melangeforce = 0
contact = 1
stabilize = 1
wall_int = 1
skip_vertex = 1

xTip= [0,0];
Rtip = xTip;
QT = eye(2);
loadstress = 'y';
Tfact = 1;
%problem flags
elemType = 'T3' ;
typeCrack = 'Static' ;
stressState = 'PlaneStrain' ;
%typeProblem = 'eCrkTen' ; %choose type of problem
typeProblem = 'ISSM' ; %choose type of problem
%typeProblem = 'yTraction' ; %choose type of problem
Kpen = 1e11;

% import rifts
% srift1 = shaperead('../../Data/2013_14_cracka_open.shp');
srift2 = shaperead('./test_rift.shp');
%srift2 = shaperead('../../Data/2013_14_crackb_open.shp');
%srift = shaperead('./Data/cracka_short_2009-10.shp');
xs = srift2.X;
ys = srift2.Y;
%xs2 = srift2.X;
%ys2 = srift2.Y;
xs(end)=[];ys(end)=[];
xs(end)=[];ys(end)=[];
xs(end)=[];ys(end)=[];
xs(1)=[];ys(1)=[];
xs(1)=[];ys(1)=[];
xs(1)=[];ys(1)=[];
xs(1)=[];ys(1)=[];

xCr(1).coor = [xs',ys'] 
%xCr(1).coor = [xs(1),ys(1);xs(4),ys(4);xs(7),ys(7)] 
%{keyboard %}


%%geometry and mesh generation
%L = 1; D = 1 ;
%rd = 0.0 ;
%ndiv(1) =  5;
%ndiv(2) = 5;
%[node,element,bcNodes,edgNodes] = createmesh(ndiv,rd) ;

load import_issm_holly1 
%element = element(1:5,:);
% we are going to use triangulation to create a Fintx...
% this is necessary if we plan on doing some refinement in the vicinity of the rift(s)
TR = triangulation(element,node);
cpos = TR.incenter;

FintX = scatteredInterpolant(cpos(:,1),cpos(:,2),ISSM_xx);
FintY = scatteredInterpolant(cpos(:,1),cpos(:,2),ISSM_yy);
FintXY = scatteredInterpolant(cpos(:,1),cpos(:,2),ISSM_xy);
FintH = scatteredInterpolant(cpos(:,1),cpos(:,2),ISSM_H');

clear global ISSM_xx ISSM_yy ISSM_xy ISSM_H% without the global these only clear in this workspace!!
%figure(1)
%triplot(TR);

% WARNING: This takes time
%plotMesh(node,element,elemType,'b-','yes',figure('visible','off'))
%print([results_path,'/original_mesh'],'-dpng','-r200')

if Hidden
  f = figure('visible','off')
else
  figure()
end

indx = find(cpos(:,1)>-2.5e5);
indy = find(cpos(:,2)<-1e6);
in = intersect(indx,indy);
size(in)
element = element(in,:);
cpos    = cpos(in,:);
f = figure();
plotMesh(node,element,elemType,'b-','yes',f)
hold on
plot(node(bc_fix,1),node(bc_fix,2),'r*')
plot(node(bc_front,1),node(bc_front,2),'cs')
print([results_path,'/section_mesh'],'-dpng','-r200')
%keyboard
clf(f);
% resize bc_fix and edge_nodes to 'in' nodes
all_node_num = unique(element);
numnode = size(all_node_num,1); 
node = node(all_node_num,:);
for i=1:3*length(element)
  element(i) = find(all_node_num==element(i));
end
edges_front2 = [];
for i=1:length(edges_front)
  if find(all_node_num==element(i,1)) & find(all_node_num==element(i,2));
    edges_front2 = [edges_front2 ; find(all_node_num==edges_front(i)),find(all_node_num==edges_front(i,2))]
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


% refinement using ameshref
path(path,'/home/antarctica/Softs/ameshref/refinement/')
in = f_find_points_xCr(cpos,xCr,35000)
%indx = find(cpos(:,1)>-20e3 & cpos(:,1)<110e3 );
%indy = find(cpos(:,2)>-1180e3 & cpos(:,2)<-1080e3 );

[node,element] = TrefineRG(node,element,in);

TR = triangulation(element,node);
triplot(TR);
print([results_path,'/mesh_refinement1'],'-dpng','-r200')
%keyboard
clf(f)

cpos = TR.incenter;
in = f_find_points_xCr(cpos,xCr,25000)


%indx = find(cpos(:,1)>0e3 & cpos(:,1)<90e3 );
%indy = find(cpos(:,2)>-1160e3 & cpos(:,2)<-1100e3 );

%in = intersect(indx,indy);
[node,element] = TrefineRG(node,element,in);

TR = triangulation(element,node);
triplot(TR);
print([results_path,'/mesh_refinement2'],'-dpng','-r200')
%keyboard
clf(f)

cpos = TR.incenter;
in = f_find_points_xCr(cpos,xCr,10000,20000)


%indx = find(cpos(:,1)>10e3 & cpos(:,1)<80e3 );
%indy = find(cpos(:,2)>-1150e3 & cpos(:,2)<-1110e3 );

%in = intersect(indx,indy);
[node,element] = TrefineRG(node,element,in);

TR = triangulation(element,node);
triplot(TR);
xl = xlim;
yl = ylim;
hold on
plot(node(bc_fix,1),node(bc_fix,2),'r*')
plot(node(bc_front,1),node(bc_front,2),'cs')
xlim(xl);
ylim(yl);
print([results_path,'/mesh_refinement2_with_bc'],'-dpng','-r200')
%keyboard
clf(f);

cpos = TR.incenter;
in = f_find_points_xCr(cpos,xCr,5000,10000)
[node,element] = TrefineRG(node,element,in);

TR = triangulation(element,node);
cpos = TR.incenter;
in = f_find_points_xCr(cpos,xCr,3000,8000)
[node,element] = TrefineRG(node,element,in);

TR = triangulation(element,node);
cpos = TR.incenter;
in = f_find_points_xCr(cpos,xCr,1000,4500)
[node,element] = TrefineRG(node,element,in);

xlim([min(xs)-30000,max(xs)+30000])
ylim([min(ys)-30000,max(ys)+30000])
print([results_path,'/mesh_refinement4'],'-dpng','-r200')
clf(f)
TR = triangulation(element,node);
triplot(TR);
cpos = TR.incenter;
figure(f)
patch('faces',element,'vertices',node,'facevertexcdata',FintY(cpos)); shading flat;
print([results_path,'/StressY_with_refinement'],'-dpng','-r200')
clf(f)

%in = f_find_points_xCr(cpos,xCr,4000,21000)
%indx = find(cpos(:,1)>39e3 & cpos(:,1)<65e3 );
%indy = find(cpos(:,2)>-1137e3 & cpos(:,2)<-1116e3 );

%in = intersect(indx,indy);
%[node,element] = TrefineRG(node,element,in);

%all_node_num = unique(element);
numnode = size(node,1); 
%node = node(all_node_num,:);
%for i=1:3*length(element)
  %element(i) = find(all_node_num==element(i));
%end
%for i=1:2*length(edges_front)
  %edges_front(i) = find(all_node_num==edges_front(i));
%end
%for i = 1:length(bc_fix)
  %edges_front(i) = find(all_node_num==bc_fix(i));
%end


numelem = size(element,1) 
TR = triangulation(element,node);
triplot(TR);
xl = [min(xs)-20000,max(xs)+20000];
yl = [min(ys)-20000,max(ys)+20000];
xlim(xl);
ylim(yl);
print([results_path,'/mesh_final'],'-dpng','-r200')
clf(f);
cpos = TR.incenter;
%if Hidden 
  %figure('visible','off')
%else
  %figure()
%end
%patch('faces',element,'vertices',node,'facevertexcdata',FintX(cpos)); shading flat;
%plotMesh(node,element,elemType,'b-','yes',figure())
%keyboard

bcNodes=[{bc_front} {1} {1} {bc_fix}]
edgNodes=[{edges_front} {1} {1} {1}]

%% Material properties and crack dimensions
E = 9.6e9; nu = 0.3; P = 1 ;
sigmato = P ;
if( strcmp(stressState,'PlaneStress') )
    C = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
else
    C = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
    Cm1 = E*0.1/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
end

x = [ -2,-0.3];
y = [-400000,-400000];

%%crack definition
deltaInc = 1000; numstep = 2;
%xCr(2).coor = [xs2',ys2'] 
xCr_orig = xCr;
typeProblem
plotmesh = 'YES' ; plotNode = 'no' ;
if Hidden 
  fmesh = figure('visible','off')
else
  fmesh = figure()
end

if( strcmp(plotmesh,'YES') )
    plotMesh(node,element,elemType,'b-',plotNode,fmesh)
    
    %crack plot
    for k=1:size(xCr,2)
        for kj = 1:size(xCr(k).coor,1)-1
            cr = plot(xCr(k).coor(kj:kj+1,1),xCr(k).coor(kj:kj+1,2),'r-') ;
            set(cr,'LineWidth',3);
        end
        for kj = 1:size(xCr(k).coor,1)
            %plot(xCr(k).coor(kj,1),xCr(k).coor(kj,2),'ro',...
                %'MarkerFaceColor',[.49 1 .63],'MarkerSize',5);
        end
    end
end
xlim(xl);
ylim(yl);
print([results_path,'/mesh_crack_start'],'-dpng','-r200')

TR = triangulation(element,node);
cpos = TR.incenter;
figure(f)
clf(f)
patch('faces',element,'vertices',node,'facevertexcdata',FintY(cpos)); shading flat;
print([results_path,'/StressYY_import'],'-dpng','-r200')
xlim(xl);
ylim(yl);
print([results_path,'/StressYY_import_zoom'],'-dpng','-r200')
clf(f)

figure(f)
patch('faces',element,'vertices',node,'facevertexcdata',FintX(cpos)); shading flat;
print([results_path,'/StressXX_import'],'-dpng','-r200')
xlim(xl);
ylim(yl);
print([results_path,'/StressXX_import_zoom'],'-dpng','-r200')
clf(f)

figure(f)
patch('faces',element,'vertices',node,'facevertexcdata',FintH(cpos)); shading flat;
print([results_path,'/Thickness_import'],'-dpng','-r200')
xlim(xl);
ylim(yl);
print([results_path,'/Thickness_import_zoom'],'-dpng','-r199')
clf(f)

global zoom_dim
zoom_dim = [xl;yl];

[Knumerical,ThetaInc,xCr] = mainXFEM(xCr,numstep,deltaInc)

%a = 3;
%C = 1.12 - 0.231*(a/D) + 10.55*(a/D)^2 - 21.72*(a/D)^3 + 30.39*(a/D)^4 ;
%KAnalytical000 = C*P*sqrt(pi*a) 
save([results_path,'/crack.mat'],'xCr','ThetaInc','Knumerical');

t = tiledlayout(2,2,'TileSpacing','Compact');
% tile 1
nexttile
plot([1:length(Knumerical{1,1})],Knumerical{1,1})
xlabel('step')
title('SIFs end 1')
legend({'K1','K2'})

nexttile
plot([1:length(Knumerical{1,2})],Knumerical{1,2})
xlabel('step')
title('SIFs end 2')
legend({'K1','K2'})

nexttile
plot([1:length(ThetaInc{1,1})],ThetaInc{1,1})
title('propagation angle, end 1')
xlabel('step')

nexttile
plot([1:length(ThetaInc{1,2})],ThetaInc{1,2})
title('propagation angle, end 2')
xlabel('step')
%plotMesh(node+dfa*[uxAna uyAna],element,elemType,'r-',plotNode)

figure_name = ['Knum_results'];
print([results_path,'/',figure_name],'-dpng','-r300')
