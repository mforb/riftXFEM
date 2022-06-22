function [Kglobal] = KmatSTAB(E_pen,enr_node,crack_node,elem_crk,type_elem,xTip,xVertex,...
    split_elem,tip_elem,vertex_elem,corner_elem,tan_elem,pos,xCrk,Kglobal,u,melange_stab)

%declare global variables here
global node element numnode numelem elemType
global E C nu
global typeProblem typeCrack
global plotmesh plotNode
global gporder numtri
global plothelp
global orig_nn
global frictionB friction_mu
global melangeforce melange contact Cm1 wall_int
global skip_branch skip_vertex
global output_file modpen stab_mu

if isempty(stab_mu)
  stab_mu = 0.1;
end
mu = friction_mu;

% we are adding the gradient dt/du to K
if strcmp(elemType,'Q4')
  corner = [1 2 3 4 1] ;
  nnode = [-1 -1;1 -1;1 1;-1 1] ;
elseif strcmp(elemType, 'T3')
  corner = [1 2 3 1] ;
  nnode = [0 0;1 0;0 1] ;
end

if skip_vertex
  elems = split_elem
else
  elems = union(split_elem,vertex_elem);
end
if ~skip_branch
  elems = union(elems,tip_elem);
end


cc = 0;
cct = 0;

for kk = 1:size(xCrk,2) %what's the crack?
  for ii=1:size(elems,1)
    iel = elems(ii) ;
    cc = cc + 1;
    sctr = element(iel,:) ;
    skip = 0;
    nn = length(sctr) ;
    if ( melange | melangeforce ) 
      [flag1,width] = f_find_melange(iel,xCrk);
      gn_lim = width;
    else
      gn_lim = 0;
    end

    [A,BrI,~,~,~] = f_enrich_assembly(iel,pos,type_elem,elem_crk,enr_node);
    [ap,apg] = f_crack_wall(iel,nnode,corner,tip_elem,vertex_elem,elem_crk,xTip,xVertex,crack_node); % elem_crk in natural coordinates
    ap = f_align_lp_gc(ap,[apg(1,:),apg(end,:)],sctr);
    if skip_branch
      if ~isempty(BrI)
        continue
      end
    end
    [W,Q] = quadrature(wall_int,'GAUSS',1) ;
    

    for seg = 1:length(ap)-1
      p = ap(seg:seg+1,:);
      pg = [apg(seg,:),apg(seg+1,:)];
      % find the distance between the two intersects (should be able to do this with det(J)
      [l,nv,mv,nnt,nmt,mmt] = f_segment_dist(pg);
      JO = l/2;

      for k_in = 1:length(Q) 
        [Np,dNdxp]=lagrange_basis('L2',Q(k_in));
        gpt = Np'*p ;
        [N,dNdxi] = lagrange_basis(elemType,gpt) ;
        pint =  N' * node(sctr,:);
        Nmat = enrNmat(N,iel,type_elem,enr_node(:,kk),elem_crk,xVertex,xTip,kk,modpen);
        gn = nv*Nmat*2*u(A) + gn_lim ;
        if gn < 0
          if k_in == 1
            if seg == 1
              cct = cct + 1;
            end
          end
          % stabalization term
          Kglobal(A,A) = Kglobal(A,A) + W(k_in)*(stab_mu*(E_pen^2)/(2*E))*((Nmat'-1/3)*nnt*(Nmat-1/3))*det(JO);
        end
        elseif melange_stab
          if flag1
            Kglobal(A,A) = Kglobal(A,A) + W(k_in)*(stab_mu*(E_pen^2)/(2*E))*((Nmat'-1/3)*nnt*(Nmat-1/3))*det(JO);
          end
        end

      end
    end
  end
end
out_str = ['The ratio of elements stabilized: ',num2str(cct/cc)];
disp(out_str);
fprintf(output_file,[out_str,'\n'])
