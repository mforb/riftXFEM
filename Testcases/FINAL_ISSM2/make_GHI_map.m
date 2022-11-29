clear element node
clear TrefineRG 
global element node
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

p1 = (ISSM_xx+ISSM_yy)/2 + sqrt( (ISSM_xx - ISSM_yy).*(ISSM_xx - ISSM_yy)/4 + ISSM_xy.*ISSM_xy);
p2 = (ISSM_xx+ISSM_yy)/2 - sqrt( (ISSM_xx - ISSM_yy).*(ISSM_xx - ISSM_yy)/4 + ISSM_xy.*ISSM_xy);

g = 9.81;
rw = 1027;
ri = 917;
F_fact = -0.5*g*ri*(1-ri/rw);
ISSM_GHI = F_fact*2*ISSM_H';

%figure();
%patch('faces',element,'vertices',node,'facevertexcdata',p1);
%cm = cbrewer2('BuPu', 256);
%colormap(cm);
%axis equal;
%shading flat;

%figure();
%patch('faces',element,'vertices',node,'facevertexcdata',ISSM_GHI);
%cm = cbrewer2('BuPu', 256);
%colormap(cm);
%axis equal;
%shading flat;

f = figure();
patch('faces',element,'vertices',node,'facevertexcdata',p1 + ISSM_GHI);
cm = cbrewer2('RdBu', 256);
colormap(flipud(cm));
cax = [-3e5,3e5];
cb = colorbar();
cb.Label.String = 'Principle Stress + GHI';
cb.FontSize = 16;
ax = gca();
ax.FontSize = 16;
axis equal;
shading flat;
caxis(cax);
ylabel('Northing (km)');
xlabel('Easting (km)');
yticks(-1300000:100000:-500000);
xticks(-600000:100000:400000);

b = f_publish_fig(f,'I');
figure_name = ['GHI_map'];
print(['./',figure_name],'-dpng','-r300')
