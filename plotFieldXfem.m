function plotFieldXfem(xCrk,pos,enrich_node,u,...
    elem_crk,vertex_elem,split_elem,tip_elem,xVertex,xTip,type_elem,ipas)

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
global Hidden Print

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
    %if ~ismember(iel,tip_elem)
    for k = 1:size(xCrk,2)
        U = [U; element_disp(iel,pos(:,k),enrich_node(:,k),u,k)];
    end
    
    %choose Gauss quadrature rules for elements
    [W,Q] = gauss_rule(iel,enrich_node,elem_crk,...
        xTip,xVertex,tip_elem,split_elem,vertex_elem,xCrk) ;
    if length(W)>1 
    %keyboard 
    end
    %loop over Gauss points
    for kk = 1:size(W,1)
        
        B = [ ] ;
        Gpt = Q(kk,:) ;
        [N,dNdxi] = lagrange_basis(elemType,Gpt) ;
        JO = node(sctr,:)'*dNdxi ;
        pt = N' * node(sctr,:);
        Gpnt = N'*node(sctr,:) ;
        %if Gpnt(1)<0
          %keyboard
        %end
        
        for k = 1:size(xCrk,2)
            B = [B xfemBmat(Gpt,iel,type_elem,enrich_node(:,k),elem_crk,xVertex,k)] ;
        end
        eps_sub = B*U ;
        
        stress(iel,kk,:) = C*eps_sub ;
        strain(iel,kk,:) = eps_sub ;
        stress_pnt = [stress_pnt; pt] ;
        %strain_pnt = [strain_pnt; pt] ;
        stress_xx = [stress_xx; stress(iel,kk,1)];
        stress_yy = [stress_yy; stress(iel,kk,2)];
        stress_xy = [stress_xy; stress(iel,kk,3)];
        vonmises  = sqrt( (stress(iel,kk,1))^2 +(stress(iel,kk,2))^2 -(stress(iel,kk,1))*(stress(iel,kk,2)) + 3*(stress(iel,kk,3))^2 );
        det = (stress(iel,kk,1))*(stress(iel,kk,2))-(stress(iel,kk,3))^2;
        
        stress_val = [stress_val; vonmises] ;
        stress_val2 = [stress_val2; det ];
        strain_val = [strain_val; strain(iel,kk,2)] ;
        
        %compute analytical stress at this point
        cracklength = 100 ;
        xcTip = xCrk(1).coor(2,:) ;
        seg = xCrk(1).coor(2,:) - xCrk(1).coor(1,:) ;
        
        [anaSigma ] = exactGriffithStress(Gpnt,xcTip,seg,P,cracklength) ;
        
        anaStress_val = [anaStress_val; anaSigma(2,1)] ;
    end
    %else
    %end
end

save([results_path,'/stress',num2str(ipas),'.mat'],'stress_xx','stress_yy','stress_xy','stress_pnt',...
     stress_val','stress_val2','node','element','elemType');

if ~Print & Hidden
  return
end

% plot(stress_pnt(:,1),stress_pnt(:,2)) ;
% tri = delaunay(stress_pnt(:,1),stress_pnt(:,2)) ;
tri = delaunay(stress_pnt(:,1),stress_pnt(:,2)) ;
% triplot(tri,stress_pnt(:,1),stress_pnt(:,2)) ;
% pause
if Hidden
  figure('visible','off');
else
  figure();
end
%figure('visible','off');
clf
hold on
patch('faces',tri,'vertices',stress_pnt,'facevertexcdata',stress_val);
title('XFEM Von Mises Stress')
shading interp
colorbar
caxis([min(0,quantile(stress_val,0.1)) round(quantile(stress_val,0.995))])
%caxis([0,8])
if Print
  figure_name = ['Von_Mises_',num2str(ipas)];
  print([results_path,'/',figure_name],'-dpng','-r300')
end

if Hidden
  figure('visible','off');
else
  figure();
end
hold on
patch('faces',tri,'vertices',stress_pnt,'facevertexcdata',stress_xx);
title('Stress XX')
shading interp
colorbar
caxis([min(0,round(quantile(stress_xx,0.05))), round(quantile(stress_xx,0.995))])
if Print
  figure_name = ['StressXX_',num2str(ipas)];
  print([results_path,'/',figure_name],'-dpng','-r300')
end

if Hidden
  figure('visible','off');
else
  figure();
end
hold on
patch('faces',tri,'vertices',stress_pnt,'facevertexcdata',stress_yy);
title('Stress YY')
shading interp
colorbar
caxis([min(0,round(quantile(stress_yy,0.05))), round(quantile(stress_yy,0.995))])
if Print
  figure_name = ['StressYY_',num2str(ipas)];
  print([results_path,'/',figure_name],'-dpng','-r300')
end

if Hidden
  figure('visible','off');
else
  figure();
end
hold on
clf
[x2,y2,z2] = griddata(stress_pnt(:,1),stress_pnt(:,2),stress_yy,[0:0.1:7],[-8:0.1:8]');
hold on
pcolor(x2,y2,z2);
colorbar
caxis([min(0,round(quantile(stress_yy,0.05))), round(quantile(stress_yy,0.995))])
if Print
  figure_name = ['StressYY_b_',num2str(ipas)];
  print([results_path,'/',figure_name],'-dpng','-r300')
end


if Hidden
  figure('visible','off');
else
  figure();
end
hold on
patch('faces',tri,'vertices',stress_pnt,'facevertexcdata',stress_val2);
title('3rd Invariant of the Stress')
shading interp
colorbar
caxis([min(0,round(quantile(stress_val2,0.05))), round(quantile(stress_val2,0.995))])
if Print
  figure_name = ['Invariant3_',num2str(ipas)];
  print([results_path,'/',figure_name],'-dpng','-r300')
end



%GPcoor = stress_pnt;
%VM = stress_yy;


%minx = 0;
%maxx = 1;

%miny = -1;
%maxy = 1 ;

%fac1 = 200;
%fac2 = 2*fac1 ;
%if Hidden
  %figure('visible','off');
%else
  %figure();
%end

%%plotting new technique
%[X,Y,Z]=griddata(GPcoor(:,1),GPcoor(:,2),VM,[minx:(maxx-minx)/fac1:maxx],[miny:(maxy-miny)/fac2:maxy]','linear');
%Z=smooth3(cat(3,Z,Z),'box',[7 7 7]);
%Z=squeeze(Z(:,:,1));
%for i=1:size(Z,1)
    %for j=1:size(Z,2)
        %if isnan(Z(i,j))
            %Z(i,j)=0;
        %end
    %end
%end
%[mumble, H, grumble] = contourf(X,Y,Z,200);
%for q=1:length(H)
    %set(H(q),'LineStyle','none')
%end
%colorbar
%axis equal
%figure_name = ['VM_method2_',num2str(ipas)];
%print([results_path,'/',figure_name],'-dpng','-r300')
