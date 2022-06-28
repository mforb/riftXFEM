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
%lname = ['./MEL1_xmas_tip1_10km','/crack1.mat']; 
lname = ['./PRES_xmas_tip1_10km','/crack1.mat']; 
load(lname)
if isfield(xCrk,'melange')
  xCrk.melange = [ 0; xCrk.melange ];
  xCrk.width = [ 0, xCrk.width ];
end
TR = triangulation(element,node);
[crackLips,flagP] = f_find_cracklips( u, xCrk, 1, [], typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
zoom_dim(1,:) = [min(xCrk.coor(:,1))-20000,max(xCrk.coor(:,1))+20000];
zoom_dim(2,:) = [min(xCrk.coor(:,2))-10000,max(xCrk.coor(:,2))+10000];
f = figure();
xlabel('test4')
triplot(TR);
hold on
axis equal;
xlim(zoom_dim(1,:));
ylim(zoom_dim(2,:));
f_plotCrack(crackLips,1e3,'r-','k-','c--')

