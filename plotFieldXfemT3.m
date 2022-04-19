function plotFieldXfem(xCrk,pos,enrich_node,crack_nodes,u,...
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

for iel=1:numelem
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

mstress = mean(stress,2);
vonmises  = sqrt( (mstress(:,1,1)).^2 +(mstress(:,1,2)).^2 -(mstress(:,1,1)).*(mstress(:,1,2)) + 3*(mstress(:,1,3).^2) );

%figure('visible','off');
if Hidden
  f = figure('visible','off');
else
  f = figure();
end
hold on
patch('faces',tri,'vertices',node,'facevertexcdata',vonmises);
title('Vonmises')
shading flat 
colorbar
try 
  caxis([min(0,quantile(vonmises,0.1)) round(quantile(vonmises,0.995))])
catch 
  mi = min(vonmises(:,1,1));
  ma = max(vonmises(:,1,1));
  sl = (ma-mi);
  if sl == 0
    sl = 1;
  end
  caxis([mi+0.2*sl,ma-0.2*sl])
end
figure_name = ['Vonmises_stress_',num2str(ipas)];
print([results_path,'/',figure_name],'-dpng','-r300')
if ~isempty(zoom_dim)
  xlim(zoom_dim(1,:))
  ylim(zoom_dim(2,:))
  figure_name = ['VonmisesZoom_stress_',num2str(ipas)];
  print([results_path,'/',figure_name],'-dpng','-r300')
end
clf(f)
% pause

figure(f);
hold on
patch('faces',tri,'vertices',node,'facevertexcdata',mstress(:,1,1));
title('Stress XX')
shading flat 
colorbar
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
clf(f);

figure(f);
hold on
patch('faces',tri,'vertices',node,'facevertexcdata',mstress(:,1,2));
title('Stress YY')
shading flat 
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
clf(f);
clf(f);

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

