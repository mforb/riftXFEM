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
path(path,'../../Testcases/Tri_simple')

%declare global variables here
global L D E nu P C sigmato
global elemType typeMesh typeProblem typeCrack stressState
global xCr deltaInc numstep numcrack
global plotmesh plotNode
global node element numnode numelem bcNodes edgNodes
global plothelp
global penalty
global results_path rift_wall_pressure zoom_dim
global epsilon
global fmesh plot_stresses
results_path = './Results1';
zoom_dim = [-1,501;2980,3020];
mkdir(results_path);
epsilon = 0.2 
penalty  = 0 
rift_wall_pressure = 0

%problem flags
elemType = 'T3' ;
typeCrack = 'Static' ;
stressState = 'PlaneStrain' ;
typeProblem = 'centre' ; %choose type of problem
plot_stresses = 1;


%geometry and mesh generation
gmsh_matlab_test
%element = element(1:42,:);
%element = element(1:270,:);

element = tricheck(node,element);
numnode = size(node,1) ;
numelem = size(element,1) ;

E = 9.6e9; nu = 0.33; P = 1.5*144.53e3 ;
sigmato = P ;
if( strcmp(stressState,'PlaneStress') )
    C = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
else
    C = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
end

%crack definition
deltaInc = 0; numstep = 1;
xCr(1).coor = [-0.1 3000 ;500 3000] ;
numcrack = size(xCr,2) ;
fmesh = figure();
TR = triangulation(element,node);
%plot the mesh before proceeding
plotmesh = 'YES' ; plotNode = 'no' ;
if( strcmp(plotmesh,'YES') )
    %plotMesh(node,element,elemType,'b-',plotNode,fmesh);
    triplot(TR);
    hold on
    
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
results_path = './Results1/';
mkdir(results_path);

[Knum1,ThetaInc1,xCr1] = mainXFEM(xCr,numstep,deltaInc); 
global crack_load
crack_load = -144.53e3;
close all

fmesh = figure();
plotmesh = 'YES' ; plotNode = 'no' ;
if( strcmp(plotmesh,'YES') )
    %plotMesh(node,element,elemType,'b-',plotNode,fmesh);
    triplot(TR);
    hold on
    
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
results_path = './Results2/';
mkdir(results_path);

[Knum2,ThetaInc2,xCr2] = mainXFEM(xCr,numstep,deltaInc); 
%fmesh = figure();
%xCr(1).coor = [-0.1 0.084764928;0.4 0.084764928] ;
%plotmesh = 'YES' ; plotNode = 'no' ;
%if( strcmp(plotmesh,'YES') )
    %plotMesh(node,element,elemType,'b-',plotNode,fmesh)
    
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

%results_path = './Results1/Through';
%mkdir(results_path);
%[Knum2,ThetaInc2,xCr2] = mainXFEM(xCr,numstep,deltaInc) 

%close all
%fmesh = figure();
%xCr(1).coor = [-0.1 0.08505;0.4 0.08505] ;
%plotmesh = 'YES' ; plotNode = 'no' ;
%if( strcmp(plotmesh,'YES') )
    %plotMesh(node,element,elemType,'b-',plotNode,fmesh)
    
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

%results_path = './Results1/Above';
%mkdir(results_path);
%[Knum3,ThetaInc3,xCr3] = mainXFEM(xCr,numstep,deltaInc) 

%disp(['Results for crack going above node:'] )
%disp(['K1: ',num2str(Knum1{2}(1))] )
%disp(['K2: ',num2str(Knum1{2}(2))] )
%disp(['Theta: ',num2str(ThetaInc1{2})] )
%disp(['----------------'])

%disp(['Results for crack going through node:'] )
%disp(['K1: ',num2str(Knum2{2}(1))] )
%disp(['K2: ',num2str(Knum2{2}(2))] )
%disp(['Theta: ',num2str(ThetaInc2{2})] )
%disp(['----------------'])

%disp(['Results for crack going below node:'] )
%disp(['K1: ',num2str(Knum3{2}(1))] )
%disp(['K2: ',num2str(Knum3{2}(2))] )
%disp(['Theta: ',num2str(ThetaInc3{2})] )
%disp(['----------------'])

%disp(['Note the cracks are slightly different'])
