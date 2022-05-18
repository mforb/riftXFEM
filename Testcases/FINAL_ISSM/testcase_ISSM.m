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
global zoom_dim
epsilon = 5 

OPT = 2; Hidden = true;

same_coords = 1

rift_wall_pressure = 0
melange = 0 
melangeforce = 0

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

% import rifts
% srift1 = shaperead('../../Data/2013_14_cracka_open.shp');
srift2 = shaperead('~/Work/Shapefiles/rift_2005.shp');
%srift2 = shaperead('../../Data/2013_14_crackb_open.shp');
%srift = shaperead('./Data/cracka_short_2009-10.shp');
xs = srift2.X
ys = srift2.Y
xs(end) = []; %get rid of trailin NaN
ys(end) = [];
xCr(1).coor = [xs',ys'] 
%xCr(1).coor = [xs(1),ys(1);xs(4),ys(4);xs(7),ys(7)] 
%{keyboard %}
results_path = './ISSM_xmas_tip1_10km';
mkdir(results_path);
path(path,'/home/antarctica/Softs/ameshref/refinement/')
run_mesh_prep
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
deltaInc = 2000; numstep = 5;
%xCr(2).coor = [xs2',ys2'] 
xCr(1).tip = [0,1];
xCr_orig = xCr;
 
typeProblem
plotmesh = 'YES' ; plotNode = 'no' ;
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

%a = 3;
%C = 1.12 - 0.231*(a/D) + 10.55*(a/D)^2 - 21.72*(a/D)^3 + 30.39*(a/D)^4 ;
%KAnalytical000 = C*P*sqrt(pi*a) 
save([results_path,'/crack.mat'],'xCr','ThetaInc','Knumerical');
make_knum

close all;
clf;
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
results_path = './ISSM_xmas_tip1_20km';
mkdir(results_path);
run_mesh_prep
zoom_dim(1,:) = [min(xCr.coor(:,1))-20000,max(xCr.coor(:,1))+20000];
zoom_dim(2,:) = [min(xCr.coor(:,2))-10000,max(xCr.coor(:,2))+10000];
[Knumerical,ThetaInc,xCr] = mainXFEM(xCr,numstep,deltaInc);
save([results_path,'/crack.mat'],'xCr','ThetaInc','Knumerical');
make_knum

close all;
clf;
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
results_path = './ISSM_xmas_tip1_30km';
mkdir(results_path);
run_mesh_prep
zoom_dim(1,:) = [min(xCr.coor(:,1))-20000,max(xCr.coor(:,1))+20000];
zoom_dim(2,:) = [min(xCr.coor(:,2))-10000,max(xCr.coor(:,2))+10000];
[Knumerical,ThetaInc,xCr] = mainXFEM(xCr,numstep,deltaInc);
save([results_path,'/crack.mat'],'xCr','ThetaInc','Knumerical');
make_knum

close all;
clf;
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
results_path = './ISSM_xmas_tip1_40km';
mkdir(results_path);
run_mesh_prep
zoom_dim(1,:) = [min(xCr.coor(:,1))-20000,max(xCr.coor(:,1))+20000];
zoom_dim(2,:) = [min(xCr.coor(:,2))-10000,max(xCr.coor(:,2))+10000];
[Knumerical,ThetaInc,xCr] = mainXFEM(xCr,numstep,deltaInc);
save([results_path,'/crack.mat'],'xCr','ThetaInc','Knumerical');
make_knum

xCr(1).tip = [1,0];
close all;
clf;
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
results_path = './ISSM_xmas_tip2_10km';
mkdir(results_path);
run_mesh_prep
zoom_dim(1,:) = [min(xCr.coor(:,1))-20000,max(xCr.coor(:,1))+20000];
zoom_dim(2,:) = [min(xCr.coor(:,2))-10000,max(xCr.coor(:,2))+10000];
[Knumerical,ThetaInc,xCr] = mainXFEM(xCr,numstep,deltaInc);
save([results_path,'/crack.mat'],'xCr','ThetaInc','Knumerical');
make_knum

close all;
clf;
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
results_path = './ISSM_xmas_tip2_20km';
mkdir(results_path);
run_mesh_prep
zoom_dim(1,:) = [min(xCr.coor(:,1))-20000,max(xCr.coor(:,1))+20000];
zoom_dim(2,:) = [min(xCr.coor(:,2))-10000,max(xCr.coor(:,2))+10000];
[Knumerical,ThetaInc,xCr] = mainXFEM(xCr,numstep,deltaInc);
save([results_path,'/crack.mat'],'xCr','ThetaInc','Knumerical');
make_knum
