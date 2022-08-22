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
global loadstress FintX FintY FintXY FintH
% global variables for conversion between two coordinate systems
global Rtip QT xTip Tfact
global ISSM_xx ISSM_yy ISSM_xy
global OPT Hidden epsilon melange melangeforce Cm1 xM rift_wall_pressure
global zoom_dim modpen modocean stab_mu
epsilon = 5 

OPT = 2; 
Hidden = true;

same_coords = 1

rift_wall_pressure = 1
melange = 0 
melangeforce = 0

global wall_int stabalize Kpen penalty contact skip_branch skip_vertex
wall_int = 2; % H is evaluated on a per element basis, therefore there is no reason to use more then interface guass point
stabalize = 1;
stab_mu = 1;
contact = 1;
Kpen = 1e12 ;
penalty = 1;
skip_branch = 0;
skip_vertex = 0;
modpen = 0;

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
xCr(1).coor = [fliplr(xs)',fliplr(ys)'] 
%xCr(1).coor = [xs(1),ys(1);xs(4),ys(4);xs(7),ys(7)] 
%{keyboard %}
results_path = './Test1';
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
deltaInc = 2500; numstep = 1;
%xCr(2).coor = [xs2',ys2'] 
xCr(1).tip = [1,0];
xCr_orig = xCr;

count = 1

for mpen = 1:3 
  for Kpen_i = 1:2 
    for stab_i = 1:3 
      for wall_ints = 1:2
 
      %wall_int = wall_ints*wall_ints; 
      wall_int = wall_ints; 
      sv = [10,1,0.01];
      stab_mu = sv(stab_i);
      %stab_mu = 0.1;
      contact = 1;
      Kv = [1e11,1e15];
      Kpen = Kv(Kpen_i) ;
      penalty = 1;
      skip_branch = 0;
      skip_vertex = 0;
      if mpen == 1
        modpen = 0;
        modocean = 0;
      elseif mpen == 2
        modpen = 1;
        modocean = 0;
      else mpen == 3
        modpen = 1;
        modocean = 1;
      end
      %a = 3;
      %C = 1.12 - 0.231*(a/D) + 10.55*(a/D)^2 - 21.72*(a/D)^3 + 30.39*(a/D)^4 ;
      %KAnalytical000 = C*P*sqrt(pi*a) 

      results_path = ['./BatchTest',num2str(count)];
      mkdir(results_path);
      % Create a table with the data and variable names
      T = table(mpen, Kpen, stab_mu, wall_int,'VariableNames', { 'mod_penalty', 'penalty_stiffness','stabalizing_param','integration_points'} );
      writetable(T, [results_path,'/Conditions.txt']);
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
      [Knumerical,ThetaInc,xCrt] = mainXFEM(xCr,numstep,deltaInc);
      save([results_path,'/crack.mat'],'xCr','ThetaInc','Knumerical');
      make_knum

      close all;
      count = 1 + count;
      end
    end
  end
end


%wall_int = 2; % H is evaluated on a per element basis, therefore there is no reason to use more then interface guass point
%stabalize = 1;
%stab_mu = 0.1;
%contact = 1;
%Kpen = 1e12 ;
%penalty = 1;
%skip_branch = 0;
%skip_vertex = 0;
%modpen = 0;
%%a = 3;
%%C = 1.12 - 0.231*(a/D) + 10.55*(a/D)^2 - 21.72*(a/D)^3 + 30.39*(a/D)^4 ;
%%KAnalytical000 = C*P*sqrt(pi*a) 

%results_path = './Test2';
%mkdir(results_path);
%if Hidden
  %fmesh = figure('visible','off');
%else
  %fmesh = figure();
%end

%if( strcmp(plotmesh,'YES') )
    %plotMesh(node,element,elemType,'b-',plotNode,fmesh)
    
    %%crack plot
    %for k=1:size(xCr,2)
        %for kj = 1:size(xCr(k).coor,1)-1
            %cr = plot(xCr(k).coor(kj:kj+1,1),xCr(k).coor(kj:kj+1,2),'r-') ;
            %set(cr,'LineWidth',3);
        %end
        %for kj = 1:size(xCr(k).coor,1)
            %%plot(xCr(k).coor(kj,1),xCr(k).coor(kj,2),'ro',...
                %%'MarkerFaceColor',[.49 1 .63],'MarkerSize',5);
        %end
    %end
%end
%[Knumerical,ThetaInc,xCrt] = mainXFEM(xCr,numstep,deltaInc);
%save([results_path,'/crack.mat'],'xCr','ThetaInc','Knumerical');
%make_knum

%close all;

%wall_int = 2; % H is evaluated on a per element basis, therefore there is no reason to use more then interface guass point
%stabalize = 1;
%stab_mu = 0.1;
%contact = 1;
%Kpen = 1e13 ;
%penalty = 1;
%skip_branch = 0;
%skip_vertex = 0;
%modpen = 0;
%%a = 3;
%%C = 1.12 - 0.231*(a/D) + 10.55*(a/D)^2 - 21.72*(a/D)^3 + 30.39*(a/D)^4 ;
%%KAnalytical000 = C*P*sqrt(pi*a) 

%results_path = './Test3';
%mkdir(results_path);
%if Hidden
  %fmesh = figure('visible','off');
%else
  %fmesh = figure();
%end

%if( strcmp(plotmesh,'YES') )
    %plotMesh(node,element,elemType,'b-',plotNode,fmesh)
    
    %%crack plot
    %for k=1:size(xCr,2)
        %for kj = 1:size(xCr(k).coor,1)-1
            %cr = plot(xCr(k).coor(kj:kj+1,1),xCr(k).coor(kj:kj+1,2),'r-') ;
            %set(cr,'LineWidth',3);
        %end
        %for kj = 1:size(xCr(k).coor,1)
            %%plot(xCr(k).coor(kj,1),xCr(k).coor(kj,2),'ro',...
                %%'MarkerFaceColor',[.49 1 .63],'MarkerSize',5);
        %end
    %end
%end
%[Knumerical,ThetaInc,xCrt] = mainXFEM(xCr,numstep,deltaInc);
%save([results_path,'/crack.mat'],'xCr','ThetaInc','Knumerical');
%make_knum

%close all;

%wall_int = 2; % H is evaluated on a per element basis, therefore there is no reason to use more then interface guass point
%stabalize = 1;
%stab_mu = 0.1;
%contact = 1;
%Kpen = 1e14 ;
%penalty = 1;
%skip_branch = 0;
%skip_vertex = 0;
%modpen = 0;

%%a = 3;
%%C = 1.12 - 0.231*(a/D) + 10.55*(a/D)^2 - 21.72*(a/D)^3 + 30.39*(a/D)^4 ;
%%KAnalytical000 = C*P*sqrt(pi*a) 

%results_path = './Test4';
%mkdir(results_path);
%if Hidden
  %fmesh = figure('visible','off');
%else
  %fmesh = figure();
%end

%if( strcmp(plotmesh,'YES') )
    %plotMesh(node,element,elemType,'b-',plotNode,fmesh)
    
    %%crack plot
    %for k=1:size(xCr,2)
        %for kj = 1:size(xCr(k).coor,1)-1
            %cr = plot(xCr(k).coor(kj:kj+1,1),xCr(k).coor(kj:kj+1,2),'r-') ;
            %set(cr,'LineWidth',3);
        %end
        %for kj = 1:size(xCr(k).coor,1)
            %%plot(xCr(k).coor(kj,1),xCr(k).coor(kj,2),'ro',...
                %%'MarkerFaceColor',[.49 1 .63],'MarkerSize',5);
        %end
    %end
%end
%[Knumerical,ThetaInc,xCrt] = mainXFEM(xCr,numstep,deltaInc);
%save([results_path,'/crack.mat'],'xCr','ThetaInc','Knumerical');
%make_knum

%close all;
%wall_int = 4; % H is evaluated on a per element basis, therefore there is no reason to use more then interface guass point
%stabalize = 1;
%stab_mu = 0.1;
%contact = 1;
%Kpen = 1e12 ;
%penalty = 1;
%skip_branch = 0;
%skip_vertex = 0;
%modpen = 1;
%%a = 3;
%%C = 1.12 - 0.231*(a/D) + 10.55*(a/D)^2 - 21.72*(a/D)^3 + 30.39*(a/D)^4 ;
%%KAnalytical000 = C*P*sqrt(pi*a) 

%results_path = './Test5';
%mkdir(results_path);
%if Hidden
  %fmesh = figure('visible','off');
%else
  %fmesh = figure();
%end

%if( strcmp(plotmesh,'YES') )
    %plotMesh(node,element,elemType,'b-',plotNode,fmesh)
    
    %%crack plot
    %for k=1:size(xCr,2)
        %for kj = 1:size(xCr(k).coor,1)-1
            %cr = plot(xCr(k).coor(kj:kj+1,1),xCr(k).coor(kj:kj+1,2),'r-') ;
            %set(cr,'LineWidth',3);
        %end
        %for kj = 1:size(xCr(k).coor,1)
            %%plot(xCr(k).coor(kj,1),xCr(k).coor(kj,2),'ro',...
                %%'MarkerFaceColor',[.49 1 .63],'MarkerSize',5);
        %end
    %end
%end
%[Knumerical,ThetaInc,xCrt] = mainXFEM(xCr,numstep,deltaInc);
%save([results_path,'/crack.mat'],'xCr','ThetaInc','Knumerical');
%make_knum
