function plotFieldXfem_pp(xCrk,pos,enrich_node,crack_nodes,u,...
    elem_crk,vertex_elem,corner_elem,split_elem,tip_elem,xVertex,xTip,type_elem,ipas)

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
            B = [B xfemBmat(Gpt,iel,type_elem,enrich_node(:,k),elem_crk,xVertex,crack_nodes,k)];
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
title('Vonmises','FontSize',fontSize1)
shading flat 
colorbar
node_vm = zeros(size(node(:,1)));
for i = 1:length(node)
  [els,~] = find(element==i);
  node_vm(i) = mean(vonmises(els));
end
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
v1  = linspace(mi,ma-0.4*sl,100);
vl = logspace(max(log10(mi),1),log10(ma-0.4*sl),20);
figure(f2);
[C,h] = tricontour(element,node(:,1),node(:,2),node_vm,v1);
figure(f3);
[C,h] = tricontour(element,node(:,1),node(:,2),node_vm,vl);

figure(f);
colorbar
title('Vonmises','FontSize',fontSize1)
caxis(ca);
figure_name = ['Vonmises_stress_',num2str(ipas)];
print([results_path,'/',figure_name],'-dpng','-r300')
figure(f2);
caxis(ca);
colorbar
colormap(cm);
title('Vonmises','FontSize',fontSize1)
figure_name = ['ContourVM_stress_lin',num2str(ipas)];
print([results_path,'/',figure_name],'-dpng','-r300');
saveas(f2,[results_path,'/',figure_name],'epsc');
figure(f3);
colorbar
caxis(ca2);
colormap(cm);
set(gca,'ColorScale','log');
title('Vonmises','FontSize',fontSize1)
figure_name = ['ContourVM_stress_log',num2str(ipas)];
print([results_path,'/',figure_name],'-dpng','-r300');
saveas(f3,[results_path,'/',figure_name],'epsc');

if ~isempty(zoom_dim)
  indx = find(cpos(:,1)>zoom_dim(1,1));
  indx = find(cpos(indx,1)<zoom_dim(1,2));
  indy = find(cpos(:,2)>zoom_dim(2,1));
  indy = find(cpos(indy,2)<zoom_dim(2,2));
  in = intersect(indx,indy);
  vm_s    = vonmises(in,1,1);
  try 
    ca = [min(0,quantile(vonmises,0.1)) round(quantile(vonmises,0.995))]
    mi = min(vm_s(:,1,1));
    ma = max(vm_s(:,1,1));
  catch 
    mi = min(vm_s(:,1,1));
    ma = max(vm_s(:,1,1));
    sl = (ma-mi);
    if sl == 0
      sl = 1;
    end
    ca = [mi,ma-0.9*sl];
    ca2 = [mi,ma];
  end
  figure(f);
  xlim(zoom_dim(1,:))
  ylim(zoom_dim(2,:))
  figure_name = ['Vonmises_stress_zoom',num2str(ipas)];
  print([results_path,'/',figure_name],'-dpng','-r300')
  figure(f2);
  xlim(zoom_dim(1,:))
  ylim(zoom_dim(2,:))
  figure_name = ['ContourVM_stress_lin_zoom',num2str(ipas)];
  print([results_path,'/',figure_name],'-dpng','-r300');
  saveas(f2,[results_path,'/',figure_name],'epsc');
  figure(f3);
  xlim(zoom_dim(1,:))
  ylim(zoom_dim(2,:))
  figure_name = ['ContourVM_stress_log_zoom',num2str(ipas)];
  print([results_path,'/',figure_name],'-dpng','-r300');
  saveas(f3,[results_path,'/',figure_name],'epsc');

  v1  = linspace(mi,ma,50);
  vl = logspace(max(log10(mi),1),log10(ma),15);
  figure(f2);

  figure(f);
  caxis(ca);
  figure_name = ['Vonmises_stress_zoom2',num2str(ipas)];
  print([results_path,'/',figure_name],'-dpng','-r300')
  figure(f2);
  clf();
  caxis(ca);
  [C,h] = tricontour(element,node(:,1),node(:,2),node_vm,v1);
  xlim(zoom_dim(1,:))
  ylim(zoom_dim(2,:))
  figure_name = ['ContourVM_stress_lin_zoom2',num2str(ipas)];
  print([results_path,'/',figure_name],'-dpng','-r300');
  saveas(f2,[results_path,'/',figure_name],'epsc');

  figure(f3);
  clf();
  [C,h] = tricontour(element,node(:,1),node(:,2),node_vm,vl);
  caxis(ca2);
  set(gca,'ColorScale','log');
  xlim(zoom_dim(1,:))
  ylim(zoom_dim(2,:))
  figure_name = ['ContourVM_stress_log_zoom2',num2str(ipas)];
  print([results_path,'/',figure_name],'-dpng','-r300');
  saveas(f3,[results_path,'/',figure_name],'epsc');

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
cm = cbrewer2('RdBu', 256);
colormap(cm);
try 
  caxis([min(0,quantile(mstress(:,1,1),0.1)) round(quantile(mstress(:,1,1),0.995))])
catch 
  mi = min(mstress(:,1,1));
  ma = max(mstress(:,1,1));
  sl = (ma-mi);
  if sl == 0
    sl = 1;
  end
  caxis([mi+0.2*sl,ma-0.2*sl])
end

figure_name = ['Stress_xx_',num2str(ipas)];
print([results_path,'/',figure_name],'-dpng','-r300')
if ~isempty(zoom_dim)
  xlim(zoom_dim(1,:))
  ylim(zoom_dim(2,:))
  figure_name = ['StressZoom_xx_',num2str(ipas)];
  print([results_path,'/',figure_name],'-dpng','-r300')
end

figure(f);
clf();
hold on
patch('faces',tri,'vertices',node,'facevertexcdata',mstress(:,1,2));
title('Stress YY','FontSize',fontSize1)
shading flat 
colorbar
cm = cbrewer2('PuOr', 256);
colormap(cm);
try 
  caxis([min(0,quantile(mstress(:,1,2),0.1)) round(quantile(mstress(:,1,2),0.995))])
catch 
  mi = min(mstress(:,1,2));
  ma = max(mstress(:,1,2));
  sl = (ma-mi);
  if sl == 0
    sl = 1;
  end
  caxis([mi+0.2*sl,ma-0.2*sl])
  warning('quantile not available- caxis boundaries are arbitrary-ish')
end
colorbar();
figure_name = ['Stress_yy_',num2str(ipas)];
print([results_path,'/',figure_name],'-dpng','-r300')
if ~isempty(zoom_dim)
  xlim(zoom_dim(1,:))
  ylim(zoom_dim(2,:))
  figure_name = ['StressZoom_yy_',num2str(ipas)];
  print([results_path,'/',figure_name],'-dpng','-r300')
end
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

