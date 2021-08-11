%to solve plate with an edge crack. the code runs for
%different mesh sizes and computes the SIF for them
%If using these codes for research or industrial purposes, please cite:
%Title: An extended finite element library
%Author(s): BORDAS, S; NGUYEN, PV; DUNANT, C; et al.
%Source: INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN ENGINEERING
%Volume: 71 Issue: 6 Pages: 703-732 Year: AUG 6 2007 DOI: 10.1002/nme.1966
%
%Title: Numerical integration over arbitrary polygonal domains based on
%Schwarz-Christofel conformal maping.
%Author(s): S Natarajan; S Bordas; D.R. Mahapatra
%Souce: INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN ENGINEERING
%DOI: 10.1002/nme.2589
%
%Author(s): Sundar, Stephane
%Other contributor(s): Stefano, Phu, David
%Modified: Jan 2011
%--------------------------------------------------------------------
clear all
close all
clc

tic

format long

%---- Define path for subroutines
path(path,'./Crackprocessing')
path(path,'./Mesh')
path(path,'./Routines_XFEM')

%declare global variables here
global L D E nu P C sigmato
global elemType typeMesh typeProblem typeCrack stressState
global xCr deltaInc numstep numcrack
global plotmesh plotNode
global node element numnode numelem bcNodes edgNodes

%problem flags
elemType = 'Q4' ;
typeCrack = 'Static' ;
stressState = 'PlaneStrain' ;
typeProblem = 'eCrkTen' ; %choose type of problem

%geometry and mesh generation
L = 2; D = 1 ;
rd = 0.0 ;
ndiv(1) = 5 ;
ndiv(2) = 2*ndiv(1) ;
[node,element,bcNodes,edgNodes] = createmesh(ndiv,rd) ;
numnode = size(node,1) ;
numelem = size(element,1) ;

E = 1e3; nu = 0.3; P = 1 ;
sigmato = P ;
if( strcmp(stressState,'PlaneStress') )
    C = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
else
    C = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
end

%crack definition
deltaInc = 0; numstep = 1;
xCr(1).coor = [-0.1 0.0;0.3 0.0] ;
numcrack = size(xCr,2) ;

%plot the mesh before proceeding
plotmesh = 'YES' ; plotNode = 'no' ;
if( strcmp(plotmesh,'YES') )
    plotMesh(node,element,elemType,'b-',plotNode)
    
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

a = 0.3 ;
C = 1.12 - 0.231*(a/D) + 10.55*(a/D)^2 - 21.72*(a/D)^3 + 30.39*(a/D)^4 ;
KAnalytical = C*P*sqrt(pi*a) 

Knumerical/KAnalytical