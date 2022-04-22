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

%define (and make) a path for results
results_path = './Tri_pressure_ocean4';
mkdir(results_path);
%copyfile('../Testcase.m',results_path);

%declare global variables here
global L D E nu P C sigmato
global Cm1 Cm2
global elemType typeMesh typeProblem typeCrack stressState
global xCr deltaInc numstep numcrack
global plotmesh plotNode
global node element numnode numelem bcNodes edgNodes
global plothelp
global penalty fixedF contact Kpen
global frictionB friction_mu epsilon FintH
global OPT Hidden epsilon melange melangeforce Cm1 xM rift_wall_pressure
global fmesh results_path same_coords
epsilon = 1e-7
same_coords = 1
plothelp = 0
rift_wall_pressure = 1
Kpen = 1e11
%problem flags
elemType = 'T3' ;
typeCrack = 'Static' ;
stressState = 'PlaneStrain' ;
typeProblem = 'eCrkTen' ; %choose type of problem
contact = 1;



%geometry and mesh generation
read_gmesh
TR = triangulation(element,node);
cpos = TR.incenter;

FintH = scatteredInterpolant(cpos(:,1),cpos(:,2),3*ones(length(cpos),1));

element = tricheck(node,element);
numnode = size(node,1) ;
numelem = size(element,1) ;

E = 1e6; nu = 0.3; P = 1 ;
sigmato = 500. ;
if( strcmp(stressState,'PlaneStress') )
    C = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
    Cm1 = E/10*(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
else
    C = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
end

%crack definition
deltaInc = 0.05; numstep = 5;
xCr(1).coor = [-0.1 -0.1;0.1, 0.1] ;
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

[Knum1,ThetaInc1,xCr1] = mainXFEM(xCr,numstep,deltaInc); 

%close all
%fixedF = [0,1];



%xCr(1).coor = [-0.1 0.084764928;0.4 0.084764928] ;
%plotmesh = 'YES' ; plotNode = 'no' ;
%if( strcmp(plotmesh,'YES') )
    %plotMesh(node,element,elemType,'b-',plotNode)
    
    %%crack plot
    %for k=1:size(xCr,2)
        %for kj = 1:size(xCr(k).coor,1)-1
            %cr = plot(xCr(k).coor(kj:kj+1,1),xCr(k).coor(kj:kj+1,2),'r-') ;
            %set(cr,'LineWidth',3);
        %end
        %for kj = 1:size(xCr(k).coor,1)
            %plot(xCr(k).coor(kj,1),xCr(k).coor(kj,2),'ro',...
                %'MarkerFaceColor',[.49 1 .63],'MarkerSize',5);
        %end
    %end
%end

%[Knum2,ThetaInc2,xCr2] = mainXFEM(xCr,numstep,deltaInc) 

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
