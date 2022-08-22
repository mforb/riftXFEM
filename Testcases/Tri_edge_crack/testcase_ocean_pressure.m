%path and variable initialisation
clear all
close all
clc
tic

format long

%---- Define path for subroutines
path(path,'../../')
path(path,'../../Crackprocessing')
path(path,'../../Mesh')
path(path,'../../Routines_XFEM')
path(path,'../../Routines_ICEM')
%------------------
results_path = './Ice_shelf_ocean_pressure';
mkdir(results_path);

%declare global variables here
global L D E nu P C sigmato
global Cm1 Cm2
global elemType typeMesh typeProblem typeCrack stressState
global xCr deltaInc numstep numcrack
global plotmesh plotNode
global node element numnode numelem bcNodes edgNodes
global plothelp
% global variables for stress
global loadstress FintX FintY FintXY FintH
% global variables for conversion between two coordinate systems
global Rtip QT xTip Tfact
global ISSM_xx ISSM_yy ISSM_xy
global OPT Hidden epsilon
global results_path rift_wall_pressure
global same_coords
global fmesh

same_coords = 1

rift_wall_pressure = 0

epsilon = 2 
plothelp = 0
contact = 0
Kpen = 1e7
penalty = 0;
%problem flags
elemType = 'T3' ;
typeCrack = 'Static' ;
stressState = 'PlaneStress' ;

OPT = 2; Hidden = true;

xTip= [0,0];
Rtip = xTip;
QT = eye(2);
loadstress = 'n';
Tfact = 1;
%problem flags
elemType = 'T3' ;
typeCrack = 'Static' ;
stressState = 'PlaneStrain' ;
typeProblem = 'eCrkTen2' ; %choose type of problem
%typeProblem = 'Test' ; %choose type of problem
%typeProblem = 'yTraction' ; %choose type of problem

% read the mesh
read_gmesh

TR = triangulation(element,node);
cpos = TR.incenter;

FintX = scatteredInterpolant(cpos(:,1),cpos(:,2),1e4*ones(size(cpos,1),1));
FintY = scatteredInterpolant(cpos(:,1),cpos(:,2),2.9e6*ones(size(cpos,1),1));
FintXY = scatteredInterpolant(cpos(:,1),cpos(:,2),0.3e3*ones(size(cpos,1),1));
FintH = scatteredInterpolant(cpos(:,1),cpos(:,2),100*ones(size(cpos,1),1));



%geometry and mesh generation
read_gmesh
%element = element(1:42,:);
%element = element(1:270,:);

element = tricheck(node,element);
numnode = size(node,1) ;
numelem = size(element,1) ;

rw = 1027;
ri = 917;
g = 9.81;

E =1e10; nu = 0.33; P = 1 ;
%sigmato = 300*g*ri*(rw-ri)/rw;
%sigmato = 2e5;
sigmato = 0.5*100*100*g*ri*(1-ri/rw);
if( strcmp(stressState,'PlaneStress') )
    C = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
    Cm1 = E/10*(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
else
    C = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
end

%crack definition
deltaInc = 100; numstep = 1;
%xCr(1).coor = [-1501 0; -500 0 ] ;
xCr(1).coor = [500 0; 1501 0 ] ;
%xCr(1).coor = [-0.2 0;0.2 0] ;
numcrack = size(xCr,2) ;
fixedF = [];

%plot the mesh before proceeding
plotmesh = 'YES' ; plotNode = 'no' ;
fmesh = figure();
if( strcmp(plotmesh,'YES') )
    plotMesh(node,element,elemType,'b-',plotNode,fmesh)
    
    %crack plot
    for k=1:size(xCr,2)
        for kj = 1:size(xCr(k).coor,1)-1
            cr = plot(xCr(k).coor(kj:kj+1,1),xCr(k).coor(kj:kj+1,2),'r-') ;
            set(cr,'LineWidth',3);
        end
        for kj = 1:size(xCr(k).coor,1)
            plot(xCr(k).coor(kj,1),xCr(k).coor(kj,2),'ro',...
                'MarkerFaceColor',[.49 1 .63],'MarkerSize',5);
        end
    end
end

[Knumerical,ThetaInc,xCr1] = mainXFEM(xCr,numstep,deltaInc); 

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
%close all
%fixedF = [0,1];



%disp(['Results for no applied force:'] )
%disp(['K1: ',num2str(Knum1(1))] )
%disp(['K2: ',num2str(Knum1(2))] )
%disp(['Theta: ',num2str(ThetaInc1)] )
%disp(['----------------'])

%disp(['Results for applied force:'] )
%disp(['K1: ',num2str(Knum2(1))] )
%disp(['K2: ',num2str(Knum2(2))] )
%disp(['Theta: ',num2str(ThetaInc2)] )

%disp(['Note the domain integrals should be modified to include crack surface stress?'])
