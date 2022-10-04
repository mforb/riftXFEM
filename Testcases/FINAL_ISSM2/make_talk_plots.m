clear all;
path(path,'../../')
path(path,'../../Crackprocessing')
path(path,'../../Mesh')
path(path,'../../Routines_XFEM')
path(path,'../../Routines_ICEM')
path(path,genpath('~/Softs/MATLAB/TOOLS/'));
fontSize1 = 14; 
fontSize2 = 12; 
mag       = 2000;

ld = dir('../FINAL_ISSM/ISSM2_xmas_tip*');
pre = '../FINAL_ISSM/';
results_path = './Talk';
mkdir(results_path);
global results_path
global zoom_dim
global Hidden
global fontSize2 fontSize1
global elemType 
Hidden = 0;
global E C nu P
global melange
global wall_int
wall_int = 1
melange = 0
melangeforce = 0
E = 9.6e9; nu = 0.3; P = 1 ;
elemType = 'T3';
sigmato = P ;
stressState = 'PlaneStrain' ;
if( strcmp(stressState,'PlaneStress') )
    C = E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
else
    C = E/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
    Cm1 = E*0.1/(1+nu)/(1-2*nu)*[1-nu nu 0; nu 1-nu 0; 0 0 (1/2)-nu];
end

% the original crack geometry
srift2 = shaperead('~/Work/Shapefiles/rift_2005.shp');
xs = srift2.X;
ys = srift2.Y;
xs(end) = []; %get rid of trailin NaN
ys(end) = [];
xCr_original.coor = [fliplr(xs)',fliplr(ys)'];

tip1 = [ ones(1,8), zeros(1,6), ones(1,2)];
tip2 = [ zeros(1,8), ones(1,6), ones(1,2)];
f = figure(2);
clf
f.Position = [ 0, 0, 1200, 800 ];
lb = -2
ub = 32
yticks(0:5:30);
% tile 1
hold on
grid on
plot([9,9],[lb,ub],'color',[30,30,30,200]/255,'linewidth',1)
plot([15,15],[lb,ub],'color',[30,30,30,200]/255,'linewidth',1)
ylim([lb,ub]);
xlim([0,16]);
xlabel('Step')
ylabel(['SIF ($\frac{MPa}{\sqrt{m}}$)'],'interpreter','latex','FontSize',14)
ax = gca();
ax.FontSize = 16;
ax.LineWidth = 1.2;
ax.Color = 'none';
ax.TickDir = 'out';
ax.TickLength = [ 0.005 0.01 ];

f = figure(3);
clf
f.Position = [ 0, 0, 1200, 800 ];
lb = -2
ub = 34
yticks(0:5:30);
% tile 1
hold on
grid on
plot([9,9],[lb,ub],'color',[30,30,30,200]/255,'linewidth',1)
plot([15,15],[lb,ub],'color',[30,30,30,200]/255,'linewidth',1)
ylim([lb,ub]);
xlim([0,16]);
xlabel('Step')
ylabel(['SIF ($\frac{MPa}{\sqrt{m}}$)'],'interpreter','latex','FontSize',14)
ax = gca();
ax.FontSize = 16;
ax.LineWidth = 1.2;
ax.Color = 'none';
ax.TickDir = 'out';
ax.TickLength = [ 0.005 0.01 ];

c1 = cbrewer2('set2',4);
c2 = cbrewer2('dark2',4);

kl1 = [];
kl2 = [];
kr1 = [];
kr2 = [];
stp = [];
cnt = 0;
% read all of the SIF values
for i = 1:length(ld)
  dname = [pre,ld(i).name];
  strs = [dname,'/crack*.mat']; 
  lns = dir(strs);
  for i = 2:length(lns)
    cnt = cnt+1;
    cmat = lns(i).name
    lname = [dname,'/',cmat]; 
    load(lname)
    kl1 = [kl1, Knum{1}(1,end)];
    kr1 = [kr1, Knum{2}(1,end)];
    kl2 = [kl2, Knum{1}(2,end)];
    kr2 = [kr2, Knum{2}(2,end)];
    stp = [stp,cnt]
    f = figure(1);
    clf
    f.Position = [ 0, 0, 1200, 700];
    hold on
    % the crack is saved after a propagation step, so we need to modify the crack to plot 
    [crackLips,flagP,elemGap] = f_find_cracklips( u, xCrk, 1, [], typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
    cext = plot(xCrk.coor(:,1),xCrk.coor(:,2))
    xCr = xCrk;
    if xCrk.tip(1)
      xCr.coor(1,:) = [];
    end
    if xCrk.tip(2)
      xCr.coor(end,:) = [];
    end
    csol = plot(xCr.coor(:,1),xCr.coor(:,2))
    set(cext,'color',[171, 22, 45, 150]/255)
    set(csol,'color',[171, 22, 45]/255)
    set([cext,csol],'LineWidth',3)
    axis equal;
    xlim([0 9e4])
    ylim([-1.15,-1.12]*1e6)
    yticks(-1150000:10000:-1120000);
    ylabel('Northing (km)');
    xlabel('Easting (km)');
    f_publish_fig(f,'r');
    figure_name = ['crack',num2str(cnt)];
    print([results_path,'/',figure_name],'-dpng')
    figure(2)
    figure_name = ['sifL',num2str(cnt)];
    pl1 = plot(stp,kl1/1e6,'color',c1(1,:),'linewidth',3)
    pl2 = plot(stp,kl2/1e6,'color',c1(2,:),'linewidth',3)
    if cnt == 1
      set([pl1,pl2],'Marker','o');
    end
    print([results_path,'/',figure_name],'-dpng')
    delete([pl1,pl2])
    figure(3)
    figure_name = ['sifR',num2str(cnt)];
    pl1 = plot(stp,kr1/1e6,'color',c2(1,:),'linewidth',3)
    pl2 = plot(stp,kr2/1e6,'color',c2(2,:),'linewidth',3)
    if cnt == 1
      set([pl1,pl2],'Marker','o');
    end
    print([results_path,'/',figure_name],'-dpng')
    delete([pl1,pl2])
  end
end

