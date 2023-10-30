%--------------------------------------------------------------------
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

%declare global variables here
global L D E nu P C sigmato
global elemType typeMesh typeProblem typeCrack stressState
global xCr deltaInc numstep numcrack
global plotmesh plotNode
global node element numnode numelem bcNodes edgNodes
global plothelp
%global output_file
global results_path
global fmesh Hidden

Hidden = false;

results_path = './tri_1';
mkdir(results_path);

plothelp = 0

%problem flags
elemType = 'T3' ;
typeCrack = 'Static' ;
stressState = 'PlaneStress' ;
typeProblem = 'eCrkTen' ; %choose type of problem

%geometry and mesh generation
read_gmesh
%keyboard
%element = [ 1 2 4 3 ]
element = tricheck(node,element);
numnode = size(node,1) ;
numelem = size(element,1) ;

E = 1e2; nu = 0; P = 1 ;
sigmato = P ;
if( strcmp(stressState,'PlaneStress') )
    C = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
else
    C = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
end

%crack definition
deltaInc = 0; numstep = 1;
xCr(1).coor = [-0.2 0.0;-0.1 0.0] ;
numcrack = size(xCr,2) ;

%plot the mesh before proceeding
if Hidden
  fmesh = figure('visible','off');
else
  fmesh = figure();
end

plotmesh = 'YES' ; plotNode = 'no' ;
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

[Knumerical,ThetaInc,xCr] = mainXFEM(xCr,numstep,deltaInc) 

%a = 0.3 ;
%C = 1.12 - 0.231*(a/D) + 10.55*(a/D)^2 - 21.72*(a/D)^3 + 30.39*(a/D)^4 ;
%KAnalytical = C*P*sqrt(pi*a) 

