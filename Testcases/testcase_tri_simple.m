%path and variable initialisation
clear all
close all
clc

tic

format long

%---- Define path for subroutines
path(path,'../')
path(path,'../Crackprocessing')
path(path,'../Mesh')
path(path,'../Routines_XFEM')
path(path,'../Testcases/Tri_simple')

%declare global variables here
global L D E nu P C sigmato
global elemType typeMesh typeProblem typeCrack stressState
global xCr deltaInc numstep numcrack
global plotmesh plotNode
global node element numnode numelem bcNodes edgNodes
global plothelp
global contact 
global results_path rift_wall_pressure
global epsilon
global melange melangeforce Cm1 xM
results_path = './Tri_simple_results';
mkdir(results_path);
epsilon = 1e-6
plothelp = 0
penalty  = 0 
rift_wall_pressure = 'n'

%problem flags
elemType = 'T3' ;
typeCrack = 'Static' ;
stressState = 'PlaneStress' ;
typeProblem = 'eCrkTen' ; %choose type of problem

melange = 0 
melangeforce = 1 
contact = 0


%geometry and mesh generation
node = [0 0 ; 1 0 ; 2 0 ; 2  1 ; 0 1 ; 0 2 ; 1 2 ; 2 2 ; 3 0 ; 3 1; 3 2;];
element = [ 1 2 5; 5 2 7; 5 7 6; 2 4 7; 2 3 4; 7 4 8; 3 9 4 ; 4 9 10 ; 4 10 11; 4 11 8 ];
botNodes = [ 1 2 3 9 ]; 
rightNodes = [ 9 10 11 ];
topNodes = [ 11 8 7 6  ]; 
leftNodes = [ 6 5 1 ]; 

botEdges = [ 1 2; 2 3; 3 9];
rightEdges = [ 9 10; 10 11];
topEdges = [ 11 8; 8 7; 7 6 ];
leftEdges = [ 6 5; 5 1];

bcNodes = {botNodes rightNodes topNodes leftNodes};
edgNodes = {botEdges rightEdges topEdges leftEdges};

element = tricheck(node,element);
numnode = size(node,1) ;
numelem = size(element,1) ;

E = 1e3; nu = 0.3; P = 1 ;
sigmato = P ;
if( strcmp(stressState,'PlaneStress') )
    C = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
    Cm1 = E*.1/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
else
    C = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
end

%crack definition
deltaInc = 0; numstep = 1;
xCr(1).coor = [-0.1 0.6; 2.2 0.6] ;
xM.coor = [-0.1 0.6; 0.1 0.6; 0.3 0.6; 0.5 0.6 ; 2.2 0.6 ];
xM.melange = [ 0 1 1 0];
xM.width = [1 5 0.1 55];


numcrack = size(xCr,2) ;
f = figure();
%plot the mesh before proceeding
plotmesh = 'YES' ; plotNode = 'no' ;
if( strcmp(plotmesh,'YES') )
    plotMesh(node,element,elemType,'b-',plotNode,f)
    
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

