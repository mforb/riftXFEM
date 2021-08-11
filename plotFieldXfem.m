	function plotFieldXfem(xCrk,pos,enrich_node,u,...
    elem_crk,vertex_elem,split_elem,tip_elem,xVertex,xTip,type_elem)

%plot stress contour.
%
%Author(s): Sundar, Stephane, Stefano, Phu, David
%Modified: Jan 2009
%--------------------------------------------------------------------

%declare global variables
global node element numnode numelem elemType 
global E nu P
global typE stressState
global plotmesh

stress_pnt =  [ ] ;
stress_val = [ ] ;
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
    
    %choose Gauss quadrature rules for elements
    [W,Q] = gauss_rule(iel,enrich_node,elem_crk,...
        xTip,xVertex,tip_elem,split_elem,vertex_elem,xCrk) ;

    
    %loop over Gauss points
    for kk = 1:size(W,1)
        B = [ ] ;
        Gpt = Q(kk,:) ;
        [N,dNdxi] = lagrange_basis(elemType,Gpt) ;
        JO = node(sctr,:)'*dNdxi ;
        pt = N' * node(sctr,:);
        Gpnt = N'*node(sctr,:) ;
        
        for k = 1:size(xCrk,2)
            B = [B xfemBmat(Gpt,iel,type_elem,enrich_node(:,k),elem_crk,xVertex,k)] ;
        end
        eps_sub = B*U ;
        
        C = (E/1-nu^2)*[1 nu 0;nu 1 0;0 0 (1-nu)/2] ;
        
        stress(iel,kk,:) = C*eps_sub ;
        strain(iel,kk,:) = eps_sub ;
        stress_pnt = [stress_pnt; pt] ;
        strain_pnt = [strain_pnt; pt] ;
        
        stress_val = [stress_val; stress(iel,kk,2)] ;
        strain_val = [strain_val; strain(iel,kk,2)] ;
        
        %compute analytical stress at this point
        cracklength = 100 ;
        xcTip = xCrk(1).coor(2,:) ;
        seg = xCrk(1).coor(2,:) - xCrk(1).coor(1,:) ;
        
        [anaSigma ] = exactGriffithStress(Gpnt,xcTip,seg,P,cracklength) ;
        
        anaStress_val = [anaStress_val; anaSigma(2,1)] ;
    end
end

% plot(stress_pnt(:,1),stress_pnt(:,2)) ;
% tri = delaunay(stress_pnt(:,1),stress_pnt(:,2)) ;
tri = delaunay(strain_pnt(:,1),strain_pnt(:,2)) ;
% triplot(tri,stress_pnt(:,1),stress_pnt(:,2)) ;
% pause


figure
hold on
for i=1:size(tri,1)
%     fill(stress_pnt(tri(i,:),1),stress_pnt(tri(i,:),2),stress_val(tri(i,:
%     ))) ;
fill(strain_pnt(tri(i,:),1),strain_pnt(tri(i,:),2),strain_val(tri(i,:))) ;
    %     patch(stress_pnt(tri(i,:),1),stress_pnt(tri(i,:),2),stress_val(tri(i,:))) ;
end
title('XFEM Stress')
shading interp
colorbar

figure
hold on
for i=1:size(tri,1)
%     fill(stress_pnt(tri(i,:),1),stress_pnt(tri(i,:),2),stress_val(tri(i,:
%     ))) ;
fill(strain_pnt(tri(i,:),1),strain_pnt(tri(i,:),2),anaStress_val(tri(i,:))) ;
    %     patch(stress_pnt(tri(i,:),1),stress_pnt(tri(i,:),2),stress_val(tri(i,:))) ;
end
title('Analytical Stress')
shading interp
colorbar

% GPcoor = stress_pnt;
% VM = stress_val;
% 
% 
% minx = 0;
% maxx = 1;
% 
% miny = -1;
% maxy = 1 ;
% 
% fac1 = 200;
% fac2 = 2*fac1 ;
% %plotting new technique
% [X,Y,Z]=griddata(GPcoor(:,1),GPcoor(:,2),VM,[minx:(maxx-minx)/fac1:maxx],[miny:(maxy-miny)/fac2:maxy]','linear');
% Z=smooth3(cat(3,Z,Z),'box',[7 7 7]);
% Z=squeeze(Z(:,:,1));
% for i=1:size(Z,1)
%     for j=1:size(Z,2)
%         if isnan(Z(i,j))
%             Z(i,j)=0;
%         end
%     end
% end
% [mumble, H, grumble] = contourf(X,Y,Z,200);
% for q=1:length(H)
%     set(H(q),'LineStyle','none')
% end
% colorbar
% axis equal