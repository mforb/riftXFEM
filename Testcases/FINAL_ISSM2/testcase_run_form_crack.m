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
%copyfile('../Testcase.m',results_path);

%declare global variables here
global L D E nu P C sigmato
global elemType typeMesh typeProblem typeCrack stressState
global xCr xCr_orig deltaInc numstep numcrack
global plotmesh plotNode
global node element numnode numelem bcNodes edgNodes
global results_path fmesh
% global variables for stress
global loadstress FintX FintY FintXY
% global variables for conversion between two coordinate systems
global Rtip QT xTip Tfact
global ISSM_xx ISSM_yy ISSM_xy
global OPT Hidden epsilon melange melangeforce Cm1 xM rift_wall_pressure
global zoom_dim same_coords
epsilon = 5 

results_path = './ISSM_xmas_test';
mkdir(results_path);
OPT = 2; Hidden = true;

same_coords = 1

rift_wall_pressure = 0
melange = 0 
melangeforce = 0

loadstress = 'y';
%problem flags
elemType = 'T3' ;
stressState = 'PlaneStrain' ;
%typeProblem = 'eCrkTen' ; %choose type of problem
typeProblem = 'ISSM' ; %choose type of problem

%% Material properties and crack dimensions
E = 9.6e9; nu = 0.3; P = 1 ;
sigmato = P ;
if( strcmp(stressState,'PlaneStress') )
    C = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
else
    C = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
    Cm1 = E*0.1/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
end

%load('ISSM_xmas_test/crack1.mat')
%xCr = xCrk;
%numnode = length(node);

load('ISSM_xmas_tip1_10km/crack.mat')
path(path,'/home/antarctica/Softs/ameshref/refinement/')
run_mesh_prep % redoes the refinement process 
xs = xCr.coor(:,1);
ys = xCr.coor(:,2);%}

%%crack definition
deltaInc = 2500; numstep = 2;
%xCr(2).coor = [xs2',ys2'] 
xCr(1).tip = [1,0];
xCr_orig = xCr;
 
typeProblem = 'ISSM';
plotmesh = 'NO' ; plotNode = 'no' ;
if Hidden
  fmesh = figure('visible','off');
else
  fmesh = figure();
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


zoom_dim(1,:) = [min(xCr.coor(:,1))-20000,max(xCr.coor(:,1))+20000];
zoom_dim(2,:) = [min(xCr.coor(:,2))-10000,max(xCr.coor(:,2))+10000];
[Knumerical,ThetaInc,xCr] = mainXFEM(xCr,numstep,deltaInc);
