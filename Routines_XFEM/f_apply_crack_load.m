function [ Fext, elem_force ] = f_apply_crack_load( enr_node,elem_crk,type_elem,xTip,xVertex,...
    split_elem,tip_elem,vertex_elem,corner_elem,crack_node,enr_domain,elem_force,pos,xCrk,Fext,fh )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Fri Mar 18 18:39:27 NZDT 2022

%declare global variables here
global node element numnode numelem elemType
global plotmesh plotNode
global gporder numtri
global plothelp 
global rift_wall_pressure wall_int skip_branch
global modocean

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
  elem_force = zeros(2,size(element,1),wall_int*2);
end


%loop over elements
elems = union(split_elem,vertex_elem);
elems = union(elems,tip_elem);

sk_br = skip_branch;
skip_branch = 0; % the force is applied everywhere consistently (skip_branch is for penalty forces)


for kk = 1:size(xCrk,2) %what's the crack?
  for ii=1:size(elems,1)
    iel = elems(ii);
    sctr = element(iel,:) ;
    nn = length(sctr) ;
    ke = 0 ;
    [ap,apg] = f_crack_wall(iel,nnode,corner,tip_elem,vertex_elem,elem_crk,xTip,xVertex,crack_node); % elem_crk in natural coordinates
    ap = f_align_lp_gc(ap,[apg(1,:),apg(end,:)],sctr);
    for seg = 1:length(ap)-1
      p = ap(seg:seg+1,:);
      pg = [apg(seg,:),apg(seg+1,:)];
      elem_force(seg,iel,1:2:wall_int*2) = elem_force(seg,iel,1:2:wall_int*2) + 2*fh;

      [W,Q] = quadrature(wall_int,'GAUSS',1) ;
      % find the distance between the two intersects (should be able to do this with det(J)
      [l,nv,mv,nnt,nmt,mmt] = f_segment_dist(pg);
      JO = l/2;

        
      skip = 0;
      nn = length(sctr) ;
      n1 = zeros(1,nn);

      [A,~,~,~,~] = f_enrich_assembly(iel,pos,type_elem,elem_crk,enr_node);

      for k_in = 1:length(Q)
        [Np,dNdxp]=lagrange_basis('L2',Q(k_in));
        gpt = Np'*p;
        [N,dNdxi] = lagrange_basis(elemType,gpt) ;
        Nmat = enrNmat(N,iel,type_elem,enr_node(:,kk),elem_crk,xVertex,xTip,kk,modocean);
        Fext(A) = Fext(A) + 2*fh*W(k_in)*det(JO)*Nmat'*nv';
      end
    end
  end
end

skip_branch = sk_br;

