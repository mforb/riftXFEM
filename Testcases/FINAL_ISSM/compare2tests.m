path(path,'../../')
path(path,'../../Crackprocessing')
path(path,'../../Mesh')
path(path,'../../Routines_XFEM')
path(path,'../../Routines_ICEM')
path(path,genpath('~/Softs/MATLAB/TOOLS/'));
fontSize1 = 14; 
fontSize2 = 12; 
mag       = 2000;

ld = dir('ISSM_xmas_tip*');
results_path = './ISSM_tests';
mkdir(results_path);
global results_path
global zoom_dim
global Hidden
global fontSize2 fontSize1
global elemType 
Hidden = 0;
global E C nu P
E = 9.6e9; nu = 0.3; P = 1 ;
elemType = 'T3';
sigmato = P ;
stressState = 'PlaneStrain' ;
lname = ['./Test4','/crack1.mat']; 
load(lname)
TR = triangulation(element,node);
[crackLips,flagP] = f_find_cracklips( u, xCrk, 1, [], typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
zoom_dim(1,:) = [min(xCrk.coor(:,1))-20000,max(xCrk.coor(:,1))+20000];
zoom_dim(2,:) = [min(xCrk.coor(:,2))-10000,max(xCrk.coor(:,2))+10000];

t = tiledlayout(1,2,'TileSpacing','Compact');
% tile 1
nexttile
xlabel('test4')
triplot(TR);
hold on
axis equal;
xlim(zoom_dim(1,:));
ylim(zoom_dim(2,:));
f_plotCrack(crackLips,1000,'r-','k-','c--')

lname = ['./Test1','/crack1.mat']; 
load(lname)
TR = triangulation(element,node);
[crackLips,flagP] = f_find_cracklips( u, xCrk, 1, [], typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
zoom_dim(1,:) = [min(xCrk.coor(:,1))-20000,max(xCrk.coor(:,1))+20000];
zoom_dim(2,:) = [min(xCrk.coor(:,2))-10000,max(xCrk.coor(:,2))+10000];

nexttile
xlabel('test1')
triplot(TR);
hold on
axis equal;
xlim(zoom_dim(1,:));
ylim(zoom_dim(2,:));
f_plotCrack(crackLips,1000,'r-','k-','c--')


%lname = ['./Test5','/crack1.mat']; 
%load(lname)
%TR = triangulation(element,node);
%[crackLips,flagP] = f_find_cracklips( u, xCrk, 1, [], typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
%zoom_dim(1,:) = [min(xCrk.coor(:,1))-20000,max(xCrk.coor(:,1))+20000];
%zoom_dim(2,:) = [min(xCrk.coor(:,2))-10000,max(xCrk.coor(:,2))+10000];

%nexttile
%xlabel('test5')
%triplot(TR);
%hold on
%axis equal;
%xlim(zoom_dim(1,:));
%ylim(zoom_dim(2,:));
%f_plotCrack(crackLips,1000,'r-','k-','c--')

%lname = ['./Test6','/crack1.mat']; 
%load(lname)
%TR = triangulation(element,node);
%[crackLips,flagP] = f_find_cracklips( u, xCrk, 1, [], typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
%zoom_dim(1,:) = [min(xCrk.coor(:,1))-20000,max(xCrk.coor(:,1))+20000];
%zoom_dim(2,:) = [min(xCrk.coor(:,2))-10000,max(xCrk.coor(:,2))+10000];

%nexttile
%xlabel('test6')
%triplot(TR);
%hold on
%axis equal;
%xlim(zoom_dim(1,:));
%ylim(zoom_dim(2,:));
%f_plotCrack(crackLips,1000,'r-','k-','c--')

%lname = ['./Test8','/crack1.mat']; 
%load(lname)
%TR = triangulation(element,node);
%[crackLips,flagP] = f_find_cracklips( u, xCrk, 1, [], typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
%zoom_dim(1,:) = [min(xCrk.coor(:,1))-20000,max(xCrk.coor(:,1))+20000];
%zoom_dim(2,:) = [min(xCrk.coor(:,2))-10000,max(xCrk.coor(:,2))+10000];

%nexttile
%triplot(TR);
%hold on
%axis equal;
%xlim(zoom_dim(1,:));
%ylim(zoom_dim(2,:));
%xlabel('test8')
%f_plotCrack(crackLips,1000,'r-','k-','c--')

%lname = ['./Test9','/crack1.mat']; 
%load(lname)
%TR = triangulation(element,node);
%[crackLips,flagP] = f_find_cracklips( u, xCrk, 1, [], typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
%zoom_dim(1,:) = [min(xCrk.coor(:,1))-20000,max(xCrk.coor(:,1))+20000];
%zoom_dim(2,:) = [min(xCrk.coor(:,2))-10000,max(xCrk.coor(:,2))+10000];

%nexttile
%triplot(TR);
%hold on
%axis equal;
%xlim(zoom_dim(1,:));
%ylim(zoom_dim(2,:));
%xlabel('test9')
%f_plotCrack(crackLips,1000,'r-','k-','c--')
