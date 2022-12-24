clear element node
path(path,genpath('~/Softs/MATLAB/TOOLS/'));
path(path,'..');
global element node
%load ./import_holly5
load ../../FINAL_ISSM/import_issm_holly1
%element = element(1:5,:);
% we are going to use triangulation to create a Fintx...
% this is necessary if we plan on doing some refinement in the vicinity of the rift(s)
TR = triangulation(element,node);
cpos = TR.incenter;

FintX = scatteredInterpolant(cpos(:,1),cpos(:,2),ISSM_xx);
FintY = scatteredInterpolant(cpos(:,1),cpos(:,2),ISSM_yy);
FintXY = scatteredInterpolant(cpos(:,1),cpos(:,2),ISSM_xy);
FintH = scatteredInterpolant(cpos(:,1),cpos(:,2),ISSM_H');

%flow = [Vx,Vy]./sqrt(Vx.*Vx+Vy.*Vy);
flow = [Vx,Vy];
IM_xx = 2*ISSM_xx+ISSM_yy;
IM_yy = 2*ISSM_yy+ISSM_xx;

p1 = (IM_xx+IM_yy)/2 + sqrt( (IM_xx - IM_yy).*(IM_xx - IM_yy)/4 + ISSM_xy.*ISSM_xy);
p2 = (IM_xx+IM_yy)/2 - sqrt( (IM_xx - IM_yy).*(IM_xx - IM_yy)/4 + ISSM_xy.*ISSM_xy);


for i = 1:length(IM_xx);
  fl = mean(flow(element(i,:),:),1);
  alpha = atan2(fl(2),fl(1));
  QT  = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];           % for the transformation to local coordinate
  sigma = [IM_xx(i),IM_yy(i),ISSM_xy(i)];
  voit2ind    = [1 3;3 2];
  stressloc   = QT*sigma(voit2ind)*QT';
  STF1(i) = stressloc(1,1);
  STF2(i) = stressloc(2,2);
end

g = 9.81;
rw = 1023;
ri = 917;
F_fact = -0.5*g*ri*(1-ri/rw);
ISSM_GHI = F_fact*ISSM_H';

%figure();
%patch('faces',element,'vertices',node,'facevertexcdata',p2/1000);
%cm = cbrewer2('BuPu', 256);
%cax = [-50,50];
%axis equal;
%shading flat;
%colormap(cm);
%cb = colorbar();
%cb.Label.String = 'Principal Stress + GHI Loading (KPa)';
%cb.FontSize = 16;
%ax = gca();
%ax.FontSize = 16;
%axis equal;
%shading flat;
%caxis(cax);
%ylabel('Northing (km)');
%xlabel('Easting (km)');
%yticks(-1300000:100000:-500000);
%xticks(-600000:100000:400000);

%figure();
%patch('faces',element,'vertices',node,'facevertexcdata',ISSM_H');
%cm = cbrewer2('BuPu', 256);
%colormap(cm);
%cb = colorbar();
%cax = [100,400];
%caxis(cax);
%axis equal;
%shading flat;
srift2 = shaperead('../../ADVECT_100A/100a_invent.shp');
%srift2 = shaperead('../../Data/2013_14_crackb_open.shp');
%srift = shaperead('./Data/cracka_short_2009-10.shp');
xs = srift2.X
ys = srift2.Y
xs(end) = []; %get rid of trailin NaN
ys(end) = [];

coor = [];
for i = 1:length(element) 
 els = element(i,:); 
 c = mean(node(els,:),1);
 coor = [coor; c];
end 
[X,Y] = meshgrid(-600000:2500:400000,-1400000:2500:-400000);
gd1 = griddata(coor(:,1),coor(:,2),p1+ISSM_GHI,X,Y);
gd2 = griddata(coor(:,1),coor(:,2),STF1'+ISSM_GHI,X,Y);


mstruct = defaultm('stereo');
mstruct.geoid =  almanac('earth','wgs84','meters'); 
mstruct.falseeasting = 0;
mstruct.falsenorthing = 0;
mstruct.origin = [-71 0 0];
mstruct = defaultm(mstruct);
key = mstruct2geotiff(mstruct);       
key.GTModelTypeGeoKey = 1;
key.GTRasterTypeGeoKey = 2;
dx = 2500;
dy = 2500;
%key = struct( ...
    %'GTModelTypeGeoKey',[], ...
        %'GTRasterTypeGeoKey',[], ...
            %'GeographicTypeGeoKey',[]);
%key.GTModelTypeGeoKey = 2;
%key.GTRasterTypeGeoKey = 2;
%key.GeographicTypeGeoKey = 4326;
R = maprefpostings([-6e5,4e5],[-14e5,-4e5],dx,dy) 
%R = makerefmat(x_ris_moa(1)+dx/2, y_ris_moa(1) + dy/2, dx, dy);
geotiffwrite('MP_GHI.tif',gd1,R,'GeoKeyDirectoryTag',key,'TiffTags',struct('Compression',Tiff.Compression.None));
geotiffwrite('MF_GHI.tif',gd2,R,'GeoKeyDirectoryTag',key,'TiffTags',struct('Compression',Tiff.Compression.None));

%f = figure();
%patch('faces',element,'vertices',node,'facevertexcdata',(p1 + ISSM_GHI)/1000);
%cm = cbrewer2('RdBu', 256);
%colormap(flipud(cm));
%cax = [-80,80];
%cb = colorbar();
%cb.Label.String = 'Principal Stress + GHI Loading (KPa)';
%cb.FontSize = 16;
%ax = gca();
%ax.FontSize = 16;
%axis equal;
%shading flat;
%caxis(cax);
%ylabel('Northing (km)');
%xlabel('Easting (km)');
%yticks(-1300000:100000:-500000);
%xticks(-600000:100000:400000);
%b = f_publish_fig(f,'I');
%figure_name = ['GHI_map'];
%print(['./',figure_name],'-dpng','-r300')
%hold on
%delete(b)
%plr = plot(xs,ys,'-');
%set(plr,'linewidth',4);
%set(plr,'color',[0,0,0]);
%b = f_publish_fig(f,'I');
%figure_name = ['GHI_rift'];
%print(['./',figure_name],'-dpng','-r300')
%%clf()

%f = figure();
%patch('faces',element,'vertices',node,'facevertexcdata',(p1)/1000);
%cm = cbrewer2('BuPu', 256);
%colormap(cm);
%cax = [0,200];
%cb = colorbar();
%cb.Label.String = 'Principal Stress + GHI Loading (KPa)';
%cb.FontSize = 16;
%ax = gca();
%ax.FontSize = 16;
%axis equal;
%shading flat;
%caxis(cax);
%ylabel('Northing (km)');
%xlabel('Easting (km)');
%yticks(-1300000:100000:-500000);
%xticks(-600000:100000:400000);
%b = f_publish_fig(f,'I');
%figure_name = ['P1_RIS'];
%print(['./',figure_name],'-dpng','-r300')

%f = figure();
%patch('faces',element,'vertices',node,'facevertexcdata',(STF1'+ISSM_GHI)/1000);
%cm = cbrewer2('RdBu', 256);
%colormap(flipud(cm));
%cax = [-80,80];
%cb = colorbar();
%cb.Label.String = 'Principal Stress + GHI Loading (KPa)';
%cb.FontSize = 16;
%ax = gca();
%ax.FontSize = 16;
%axis equal;
%shading flat;
%caxis(cax);
%ylabel('Northing (km)');
%xlabel('Easting (km)');
%yticks(-1300000:100000:-500000);
%xticks(-600000:100000:400000);
%b = f_publish_fig(f,'I');
%figure_name = ['STF1_GHI'];
%print(['./',figure_name],'-dpng','-r300')

%f = figure();
%patch('faces',element,'vertices',node,'facevertexcdata',(STF2'+ISSM_GHI)/1000);
%cm = cbrewer2('RdBu', 256);
%colormap(flipud(cm));
%cax = [-80,80];
%cb = colorbar();
%cb.Label.String = 'Principal Stress + GHI Loading (KPa)';
%cb.FontSize = 16;
%ax = gca();
%ax.FontSize = 16;
%axis equal;
%shading flat;
%caxis(cax);
%ylabel('Northing (km)');
%xlabel('Easting (km)');
%yticks(-1300000:100000:-500000);
%xticks(-600000:100000:400000);
%b = f_publish_fig(f,'I');
%figure_name = ['STF2_GHI'];
%print(['./',figure_name],'-dpng','-r300')
