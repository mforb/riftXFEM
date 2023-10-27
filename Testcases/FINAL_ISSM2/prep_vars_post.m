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
global zoom_dim 
epsilon = 5 

OPT = 2; Hidden = true;

same_coords = 1

xTip= [0,0];
Rtip = xTip;
QT = eye(2);
loadstress = 'y';
Tfact = 1;
%problem flags
elemType = 'T3' ;
typeCrack = 'Static' ;
%stressState = 'PlaneStrain' ;
stressState = 'PlaneStress' ;
%typeProblem = 'eCrkTen' ; %choose type of problem
typeProblem = 'ISSM' ; %choose type of problem
%typeProblem = 'yTraction' ; %choose type of problem


srift2 = shaperead('./o_rift.shp');
%srift2 = shaperead('../../Data/2013_14_crackb_open.shp');
%srift = shaperead('./Data/cracka_short_2009-10.shp');
xs = srift2.X
ys = srift2.Y
xs(end) = []; %get rid of trailin NaN
ys(end) = [];
xCr(1).coor = [xs',ys'] 

clear element node
clear TrefineRG 
load ../FINAL_ISSM/import_issm_holly1
%element = element(1:5,:);
% we are going to use triangulation to create a Fintx...
% this is necessary if we plan on doing some refinement in the vicinity of the rift(s)
TR = triangulation(element,node);
cpos = TR.incenter;

FintX = scatteredInterpolant(cpos(:,1),cpos(:,2),ISSM_xx);
FintY = scatteredInterpolant(cpos(:,1),cpos(:,2),ISSM_yy);
FintXY = scatteredInterpolant(cpos(:,1),cpos(:,2),ISSM_xy);
FintH = scatteredInterpolant(cpos(:,1),cpos(:,2),ISSM_H');
clear global ISSM_xx ISSM_yy ISSM_xy % without the global these only clear in this workspace!!

