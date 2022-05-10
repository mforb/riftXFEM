function [Kglobal,Gint,elem_force] = KTmatXFEM(E_pen,enr_node,crack_node,elem_crk,type_elem,xTip,xVertex,...
    split_elem,tip_elem,vertex_elem,corner_elem,tan_elem,elem_force,pos,xCrk,xM,Kglobal,u)

%declare global variables here
global node element numnode numelem elemType
global E C nu
global typeProblem typeCrack
global plotmesh plotNode
global gporder numtri
global plothelp
global orig_nn
global frictionB friction_mu
global melangeforce contact Cm1
global wall_int
global rift_wall_pressure output_file
global skip_branch skip_vertex

mu = friction_mu;

% we are adding the gradient dt/du to K
if strcmp(elemType,'Q4')
  corner = [1 2 3 4 1] ;
  nnode = [-1 -1;1 -1;1 1;-1 1] ;
elseif strcmp(elemType, 'T3')
  corner = [1 2 3 1] ;
  nnode = [0 0;1 0;0 1] ;
end

%if skip_vertex
  %elems = union(split_elem,vertex_elem);
%else
  %elems = split_elem
%end
%if ~skip_branch
  %elems = union(elems,tip_elem);
%end
elems = split_elem
for kk = 1:size(xCrk,2) %what's the crack?
  for ii=1:size(elems,1)
    cc = cc + 1;
    iel = elems(ii) ;
    sctr = element(iel,:) ;
    skip = 0;
    nn = length(sctr) ;

    [A,BrI,QT,Tip,alpha] = f_enrich_assembly(iel,pos,type_elem,elem_crk,enr_node);
    [ap,apg,pos_l,neg_l] = f_nitsche_wall(iel,nnode,corner,tip_elem,vertex_elem,elem_crk,xTip,xVertex,crack_node); % elem_crk in natural coordinates
    ap = f_align_lp_gc(ap,apg,sctr);
    pos_g = node(sctr(pos_l),:);
    neg_g = node(sctr(neg_l),:);apg);
    if size(pos_g,1)>1
      pos_g = f_align_lp_gc(pos_g,apg,sctr);
    end
    if size(neg_g,1)>1
      neg_g = f_align_lp_gc(neg_g,apg,sctr);
    end
    pos_g = [pos_g;apg(2,:);apg(:,1);]
    neg_g = [neg_g;apg(2,:);apg(:,1);]
    measP = abs(polyarea(pos_g(:,1),pos_g(:,2)));
    measN = abs(polyarea(neg_g(:,1),neg_g(:,2)));
    
    gamma_1 = measP/(measP+measN);
    gamma_2 = measN/(measP+measN);

    sctrB = assembly(iel,enrich_node(:,kk),pos(:,kk),kK,crack_node) ;
    [A,~,~,~,~] = f_enrich_assembly(iel,pos,type_elem,elem_crk,enr_node);
     
    for seg = 1:length(ap)-1
      p = ap(seg:seg+1,:);
      pg = [apg(seg,:),apg(seg+1,:)];
      [W,Q] = quadrature(wall_int,'GAUSS',1) ;
      % find the distance between the two intersects (should be able to do this with det(J)
      [l,nv,mv,nnt,nmt,mmt] = f_segment_dist(pg);
      nx = -nv
      JO = l/2;
      alpha = 2 * l *det(C)*((gamma_1^2)/measP + (gamma_2^2)/measN);
      d2gamd2alp = 4*l*det(C)*(1/measP+1/measN);
      am = alpha*eye(2);

      for k_in = 1:length(Q)
        [N1,dNdx1]=lagrange_basis('L2',Q(k_in));
        gpt = N1'*p;
    % find the distance between the two intersects (should be able to do this with det(J)
        [N,dNdxi] = lagrange_basis(elemType,gpt) ;
        pint =  N' * node(sctr,:);
        B = xfemBmat(Gpt,iel,type_elem,enrich_node(:,k),elem_crk,xVertex,crack_node,k);
        Nmat = enrNmat(N,iel,type_elem,enr_node(:,kk),elem_crk,xVertex,kk,true);
        [NmatP,NmatN,BmatP,BmatN] = nitscheNBmat(N,pos_l,iel,type_elem,enr_node(:,kk),elem_crk,xVertex,kk);
        gn = nv*Nmat*2*u(A);
        if gn < 0
          if k_in ==1
            if seg == 1
              cct = cct + 1;
            end
          end
          Kd1 = (NmatP'*am*NmatP - gamma_1*NmatP'*nv'*C*BmatP)*W(k_in)*det(JO); 
          Kd2 = (NmatN'*am*NmatN - gamma_2*NmatN'*nx'*C*BmatN)*W(k_in)*det(JO); 
          Ko1 = (-1*NmatP'*am*NmatP - gamma_2*NmatN'*nx'*C*BmatN)*W(k_in)*det(JO); 
          Ko2 = (-1*NmatN'*am*NmatN - gamma_1*NmatP'*nv'*C*BmatP)*W(k_in)*det(JO); 
          f1n = gamma_1*C*BmatP*nv - alpha*N'
          elem_force(seg,iel,2*k_in-1) = -1*E_pen*gn;
          Gint(A) = Gint(A) + (E_pen*gn)*W(k_in)*det(JO)*Nmat'*nv';
          Kglobal(A,A) = Kglobal(A,A) + 2*E_pen*W(k_in)*Nmat'*nnt*Nmat*det(JO) ;
        end
      end
          
      if plothelp
        % this needs fixing
          Ppoint =  N' * node(sctr,:);
          Pvect = -1*det(JO)*nv;
          Np = node(sctr,:);
          Nvect = -1*det(JO)*N*nv;
          figure(2)
          %plotMesh(node,element(iel,:),'T3','r-','no')
          %if ismember(iel,split_elem)   

            %ppl = plot(Ppoint(1),Ppoint(2),'*m','linestyle','none','markersize',2)
            %%keyboard
            %delete(ppl)
          %end
          plot(Ppoint(1),Ppoint(2),'*c','linestyle','none','markersize',1)
          quiver(Ppoint(1),Ppoint(2),Pvect(1),Pvect(2),'r')
          quiver(Np(:,1),Np(:,2),Nvect(:,1),Nvect(:,2),'g')
          if frictionB
            Fvect = -1*det(JO)*mu*mv;
            NFvect = -1*det(JO)*N*mv;
            quiver(Ppoint(1),Ppoint(2),Fvect(1),Fvect(2),'k')
            quiver(Np(:,1),Np(:,2),NFvect(:,1),NFvect(:,2),'o')
          end
      end
    end
  end
end
out_str = ['The ratio of contact elements : ',num2str(cct/cc)];
disp(out_str);
fprintf(output_file,[out_str,'\n'])
