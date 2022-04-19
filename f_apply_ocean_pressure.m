function [ Fext, elem_force ] = f_apply_ocean_pressure( enrich_node,elem_crk,type_elem,xTip,xVertex,...
    split_elem,tip_elem,vertex_elem,corner_elem,crack_node,enr_domain,elem_force,pos,xCrk,Fext )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Fri Mar 18 18:39:27 NZDT 2022

%declare global variables here
global node element numnode numelem elemType
global plotmesh plotNode
global gporder numtri
global plothelp
global rift_wall_pressure

if strcmp(elemType,'Q4') 
  intType = 'GAUSS' ;
  corner = [1 2 3 4 1] ;
  nnode = [-1 -1;1 -1;1 1;-1 1] ;
else
  intType = 'TRIANGULAR';
  corner = [1 2 3 1] ;
  nnode = [0 0;1 0;0 1] ;
end
[W2,Q2] = quadrature(2,intType,2) ;
% number of non-enriched df
dfn = size(W2,1)*2;

if isempty(elem_force)
  elem_force = zeros(size(element,1),4);
end


%loop over elements
elems = union(split_elem,vertex_elem);
elems = union(elems,tip_elem);

for kk = 1:size(xCrk,2) %what's the crack?
  for ii=1:size(elems,1)
    iel = elems(ii);
    sctr = element(iel,:) ;
    nn = length(sctr) ;
    ke = 0 ;
    p = f_crack_wall(iel,nnode,corner,tip_elem,vertex_elem,elem_crk,xTip,crack_node); % elem_crk in natural coordinates
    fh = f_getHeightF(iel);
    elem_force(iel,[1,3]) = elem_force(iel,[1,3]) + fh;

    [W,Q] = quadrature(2,'GAUSS',1) ;
    [N1,dNdx1]=lagrange_basis('L2',Q(1));
    [N2,dNdx2]=lagrange_basis('L2',Q(2));
    gpts = [N1'*p; N2'*p];
    % find the distance between the two intersects (should be able to do this with det(J)
    x0 = elem_crk(iel,1) ; y0 = elem_crk(iel,2) ;
    x1 = elem_crk(iel,3) ; y1 = elem_crk(iel,4) ;
    l = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0)) ;
    nv = [(y0-y1),(x1-x0)]./l;
    mv = [(x1-x0),(y1 - y0)]./l;
    nnt = nv'*nv;
    nmt = mv'*nv; 

    JO = l/2;
      
    skip = 0;
    nn = length(sctr) ;
    n1 = zeros(1,nn);

    [A,BrI,QT,Tip,alpha] = f_enrich_assembly(iel,pos,type_elem,elem_crk,enrich_node);

    for k_in = 1:2
      gpt = gpts(k_in,:) ;
      [N,dNdxi] = lagrange_basis(elemType,gpt) ;
      pint =  N' * node(sctr,:);
      
      if any(enrich_node(sctr)==1)
        xp    = QT*(pint-Tip)';           % local coordinates
        r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
        theta = atan2(xp(2),xp(1));
        
        if ( theta > pi | theta < -pi)
            disp (['something wrong with angle ',num2str(thet)]);
        end
        if abs(abs(theta) - pi) < 0.001
         [Br_u,dBdx,dbdy] = branch_gp(r,pi,alpha);
         [Br_d,dBdx,dbdy] = branch_gp(r,-1*pi,alpha);
        else 
         [Br_u,dbdx,dbdy] = branch_gp(r,theta,alpha);
         Br_d = Br_u;
        end
      end
      nA = [1,2];
      for ni = 1:nn
        if n1(ni)
          n_row = sum(n1(1:ni));
          for i = 1:4
            Fext(A(nA)) = Fext(A(nA)) + fh*N(ni)*(Br_u(i)-BrI(n_row,i))*W(k_in)*det(JO)*nv';
            nA = [nA(1)+2,nA(2)+2];
          end
        else
          Fext(A(nA)) = Fext(A(nA)) + fh*W(k_in)*det(JO)*N(ni)*nv';
          nA = [nA(1)+2,nA(2)+2];
        end 
      end
    end
  end
end

