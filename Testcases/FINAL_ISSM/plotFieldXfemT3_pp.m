function [ca,cax,cay] = plotFieldXfem_pp(xCrk,pos,enrich_node,crack_nodes,u,...
    elem_crk,vertex_elem,corner_elem,split_elem,tip_elem,xVertex,xTip,type_elem,ipas,varargin)

%plot stress contour.
%
%Author(s): Sundar, Stephane, Stefano, Phu, David
%Modified: Jan 2009
%--------------------------------------------------------------------

%declare global variables
global node element numnode numelem elemType 
global E C nu P
global typE stressState
global plotmesh
global results_path
global ISSM_yy
global Hidden
global zoom_dim
global fontSize1 fontSize2

stress_pnt =  [ ] ;
stress_xx = [ ] ;
stress_xy = [ ] ;
stress_yy = [ ] ;
stress_val = [ ] ;
stress_val2 = [ ] ;
strain_pnt = [ ] ;
strain_val = [ ] ;
anaStress_val = [ ] ;
fac = 0 ;
ca = []; cax = []; cay = [];
if length(varargin)>0
  ca = varargin{1};
  if length(varargin)>1
    cax = varargin{2};
    if length(varargin)==3
      cay = varargin{3};
    end
  end
end


for iel=1:size(element,1)
    sctr = element(iel,:) ;
    nn = length(sctr) ;
    U = [ ];
    for k = 1:size(xCrk,2)
        U = [U; element_disp(iel,pos(:,k),enrich_node(:,k),u,k)];
    end
    if any(isnan(U))
      warning(['NaN values in enrichments of element ',num2str(iel)])
      U(isnan(U)) = 0;
    end
    %choose Gauss quadrature rules for elements
    [W,Q] = gauss_rule(iel,enrich_node,elem_crk,...
        xTip,xVertex,tip_elem,split_elem,vertex_elem,corner_elem,xCrk) ;

    for kk = 1:size(W,1)
    
    B = [ ] ;
    Gpt = Q(kk,:) ;
    [N,dNdxi] = lagrange_basis(elemType,Gpt) ;
    JO = node(sctr,:)'*dNdxi ;
    pt = N' * node(sctr,:);
    Gpnt = N'*node(sctr,:) ;
      
        for k = 1:size(xCrk,2)
            B = [B xfemBmat(Gpt,iel,type_elem,enrich_node(:,k),elem_crk,xVertex,xTip,crack_nodes,k)];
        end
        eps_sub = B*U ;
        
        stress(iel,kk,:) = C*eps_sub ;
        strain(iel,kk,:) = eps_sub ;
    end
end
tri = element;
TR = triangulation(element,node);
cpos = TR.incenter;

mstress = mean(stress,2);
vonmises  = sqrt( (mstress(:,1,1)).^2 +(mstress(:,1,2)).^2 -(mstress(:,1,1)).*(mstress(:,1,2)) + 3*(mstress(:,1,3).^2) );

%figure('visible','off');
if Hidden
  f = figure('visible','off');
  f2 = figure('visible','off');
  f3 = figure('visible','off');
else
  f = figure();
  f2 = figure();
  f3 = figure();
end
hold on
figure(f);
patch('faces',tri,'vertices',node,'facevertexcdata',vonmises);
cm = cbrewer2('BuPu', 256);
colormap(cm);
axis equal;
title('Vonmises','FontSize',fontSize1)
shading flat 
colorbar
node_vm = zeros(size(node(:,1)));
for i = 1:length(node)
  [els,~] = find(element==i);
  node_vm(i) = mean(vonmises(els));
end
if isempty(ca)
  try 
    ca = [min(0,quantile(vonmises,0.1)) round(quantile(vonmises,0.995))];
    mi = ca(1);
    ma = ca(2);
  catch 
    mi = min(vonmises(:,1,1));
    ma = max(vonmises(:,1,1));
    sl = (ma-mi);
    if sl == 0
      sl = 1;
    end
    ca = [mi,ma-0.9*sl];
    ca2 = [mi,ma];
  end
else
  mi = ca(1);
  sl = 10*(ca(2)- ca(1));
  ma = ca(1)+sl;
  ca2 = [mi,ma];
end
vl1 = logspace(3,log10(5e6),25);
vl2 = logspace(2,7,6);
figure(f2);
[C,h] = tricontf(node(:,1),node(:,2),element,node_vm,vl1);
set(h,'edgecolor',[0.1 0.1 0.1]);
set(h,'edgealpha',0.5);
axis equal;
colormap(cm);
caxis([vl1(1),vl1(end)]);
writematrix(vl1,'levels.dat','Delimiter',',')
set(gca,'ColorScale','log');
figure(f3);
[C,h] = tricontf(node(:,1),node(:,2),element,node_vm,vl2);
set(h,'edgecolor',[0.1 0.1 0.1]);
set(h,'edgealpha',0.5);
colormap(cm);
caxis([vl2(1),vl2(end)]);
axis equal;
set(gca,'ColorScale','log');

figure(f);
colorbar
title('Vonmises','FontSize',fontSize1)
caxis(ca);
figure_name = ['Vonmises_stress_',num2str(ipas)];
print([results_path,'/',figure_name],'-dpng','-r300')
caxis([0,3e4]);
figure_name = ['Vonmises_stress2_',num2str(ipas)];
print([results_path,'/',figure_name],'-dpng','-r300')

if ~isempty(zoom_dim)
  figure(f2);
  xlim(zoom_dim(1,:))
  ylim(zoom_dim(2,:))
  colorbar;
  figure_name = ['ContourVM_stress_log1_zoom',num2str(ipas)];
  print([results_path,'/',figure_name],'-dpng','-r300');
  %saveas(f2,[results_path,'/',figure_name],'epsc');
  [C,h] = tricont(node(:,1),node(:,2),element,node_vm,vl1);
  clabel(C);
  xlim([zoom_dim(1,1)-10000,zoom_dim(1,2)+10000])
  ylim([zoom_dim(2,1)-10000,zoom_dim(2,2)+10000])
  figure_name = ['ContourVM_labels1_zoom',num2str(ipas)]; 
  print([results_path,'/',figure_name],'-dpng','-r500');


  figure(f3);
  xlim(zoom_dim(1,:))
  ylim(zoom_dim(2,:))
  colorbar;
  figure_name = ['ContourVM_stress_log2_zoom',num2str(ipas)];
  print([results_path,'/',figure_name],'-dpng','-r300');
  %saveas(f3,[results_path,'/',figure_name],'epsc');
  [C,h] = tricont(node(:,1),node(:,2),element,node_vm,vl2);
  clabel(C);
  xlim([zoom_dim(1,1)-10000,zoom_dim(1,2)+10000])
  ylim([zoom_dim(2,1)-10000,zoom_dim(2,2)+10000])
  figure_name = ['ContourVM_labels2_zoom',num2str(ipas)]; 
  print([results_path,'/',figure_name],'-dpng','-r500');
end

close(f2);
close(f3);

figure(f);
clf();
hold on
patch('faces',tri,'vertices',node,'facevertexcdata',mstress(:,1,1));
title('Stress XX','FontSize',fontSize1)
shading flat 
colorbar
cm = cbrewer2('RdYlBu', 256);
colormap(cm);
%if isempty(cax)
%try 
  %cax = [min(0,quantile(mstress(:,1,1),0.1)) round(quantile(mstress(:,1,1),0.99))];
%catch 
  %mi = min(mstress(:,1,1));
  %ma = max(mstress(:,1,1));
  %sl = (ma-mi);
  %if sl == 0
    %sl = 1;
  %end
  %cax = [mi+0.2*sl,ma-0.2*sl];
%end
%end
cax = [-2e4,2e4];
caxis(cax);
axis equal;

figure_name = ['Stress_xx_',num2str(ipas)];
print([results_path,'/',figure_name],'-dpng','-r300')
%if ~isempty(zoom_dim)
  %xlim(zoom_dim(1,:))
  %ylim(zoom_dim(2,:))
  %figure_name = ['StressZoom_xx_',num2str(ipas)];
  %print([results_path,'/',figure_name],'-dpng','-r300')
%end

figure(f);
clf();
hold on
patch('faces',tri,'vertices',node,'facevertexcdata',mstress(:,1,2));
title('Stress YY','FontSize',fontSize1)
shading flat 
colorbar
cm = cbrewer2('RdYlGn', 256);
colormap(cm);
%if isempty(cay)
%try 
  %cay = [min(0,quantile(mstress(:,1,2),0.1)) round(quantile(mstress(:,1,2),0.99))];
%catch 
  %mi = min(mstress(:,1,2));
  %ma = max(mstress(:,1,2));
  %sl = (ma-mi);
  %if sl == 0
    %sl = 1;
  %end
  %cay = [mi+0.2*sl,ma-0.2*sl];
  %warning('quantile not available- caxis boundaries are arbitrary-ish')
%end
%end
cay = [-2e4,2e4];
caxis(cay);
colorbar();
axis equal;
figure_name = ['Stress_yy_',num2str(ipas)];
print([results_path,'/',figure_name],'-dpng','-r300')
%if ~isempty(zoom_dim)
  %xlim(zoom_dim(1,:))
  %ylim(zoom_dim(2,:))
  %figure_name = ['StressZoom_yy_',num2str(ipas)];
  %print([results_path,'/',figure_name],'-dpng','-r300')
%end
clf(f); close(f);

%figure;
%%figure('visible','off');
%clf
%hold on
%patch('faces',tri,'vertices',node,'facevertexcdata',ISSM_yy);
%title('Stress YY')
%shading flat 
%colorbar
%caxis([min(0,quantile(ISSM_yy,0.1)) round(quantile(ISSM_yy,0.995))])
%figure_name = ['ISSM_yy_',num2str(ipas)];
%print([results_path,'/',figure_name],'-dpng','-r300')

%figure;
%%figure('visible','off');
%clf
%hold on
%patch('faces',tri,'vertices',node,'facevertexcdata',mstress(:,1,3));
%title('Stress XY')
%shading flat 
%colorbar
%caxis([min(0,quantile(mstress(:,1,3),0.1)) round(quantile(mstress(:,1,3),0.995))])
%figure_name = ['Stress_xy_',num2str(ipas)];
%print([results_path,'/',figure_name],'-dpng','-r300')

