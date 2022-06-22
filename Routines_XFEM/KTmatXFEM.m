function [Kglobal,Gint,elem_force,ratio_c] = KTmatXFEM(E_pen,enr_node,crack_node,elem_crk,type_elem,xTip,xVertex,...
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
global melangeforce melange contact Cm1
global wall_int
global rift_wall_pressure output_file
global skip_branch skip_vertex modpen
global elem_tan

mu = friction_mu;

% we are adding the gradient dt/du to K
if strcmp(elemType,'Q4')
  corner = [1 2 3 4 1] ;
  nnode = [-1 -1;1 -1;1 1;-1 1] ;
elseif strcmp(elemType, 'T3')
  corner = [1 2 3 1] ;
  nnode = [0 0;1 0;0 1] ;
end

if plothelp
  figure(2)
  clf
  hold on
  split_nodes = find(enr_node == 2);
  tip_nodes   = find(enr_node == 1);
  n1 = plot(node(split_nodes,1),node(split_nodes,2),'r*');
  n2 = plot(node(tip_nodes,1),node(tip_nodes,2),'rs');
  set(n1,'MarkerSize',5);
  set(n2,'MarkerSize',5);
  plotMesh_numbered(node,element,elemType,'b-','no')
end
elems = union(split_elem,vertex_elem);
mel_elems = [];
for kk = 1:size(xCrk,2)
  for i=1:length(elems)                     %loop on elems (=elements selected for enrichment)
    iel = elems(i);
    [flag1,width] = f_find_melange(iel,xCrk(kk));
    if flag1
      mel_elems = [mel_elems; kk, iel, width];
    end
  end
end

% Gint
Gint =  zeros(size(u));
cc = 0;
cct = 0;

if skip_vertex
  elems = split_elem
else
  elems = union(split_elem,vertex_elem);
end
if ~skip_branch
  elems = union(elems,tip_elem);
end

if contact
  for kk = 1:size(xCrk,2) %what's the crack?
    for ii=1:size(elems,1)
      iel = elems(ii) ;
      if ( melange | melangeforce ) 
        [flag1,width] = f_find_melange(iel,xCrk);
        gn_lim = width
      else
        gn_lim = 0;
      end

      cc = cc + 1;
      sctr = element(iel,:) ;
      skip = 0;
      nn = length(sctr) ;

      [A,~,~,~,~] = f_enrich_assembly(iel,pos,type_elem,elem_crk,enr_node);
      [ap,apg] = f_crack_wall(iel,nnode,corner,tip_elem,vertex_elem,elem_crk,xTip,xVertex,crack_node); % elem_crk in natural coordinates
      ap = f_align_lp_gc(ap,[apg(1,:),apg(end,:)],sctr);
      for seg = 1:length(ap)-1
        p = ap(seg:seg+1,:);
        pg = [apg(seg,:),apg(seg+1,:)];
        [W,Q] = quadrature(wall_int,'GAUSS',1) ;
        % find the distance between the two intersects (should be able to do this with det(J)
        [l,nv,mv,nnt,nmt,mmt] = f_segment_dist(pg);
        JO = l/2;

        for k_in = 1:length(Q)
          [N1,dNdx1]=lagrange_basis('L2',Q(k_in));
          gpt = N1'*p;
      % find the distance between the two intersects (should be able to do this with det(J)
          [N,dNdxi] = lagrange_basis(elemType,gpt) ;
          pint =  N' * node(sctr,:);
          Nmat = enrNmat(N,iel,type_elem,enr_node(:,kk),elem_crk,xVertex,xTip,kk,modpen);
          gn = nv*Nmat*2*u(A) + gn_lim;
          if gn < 0 
            if k_in ==1
              if seg == 1
                cct = cct + 1;
              end
            end
            elem_force(seg,iel,2*k_in-1) = -1*E_pen*gn;
            Gint(A) = Gint(A) + (E_pen*gn)*W(k_in)*det(JO)*Nmat'*nv';
            Kglobal(A,A) = Kglobal(A,A) + 2*E_pen*W(k_in)*Nmat'*nnt*Nmat*det(JO) ;
          else          
            elem_force(seg,iel,2*k_in-1) = 0;
          end
        end
      end
    end
  end
end
ratio_c = cct/cc;

if melangeforce
  for ii=1:size(mel_elems,1)
    iel = mel_elems(ii,2) ;
    sctr = element(iel,:) ;
    nn = length(sctr) ;

    [A,~,~,~,~] = f_enrich_assembly(iel,pos,type_elem,elem_crk,enr_node);
    [ap,apg] = f_crack_wall(iel,nnode,corner,tip_elem,vertex_elem,elem_crk,xTip,xVertex,crack_node); % elem_crk in natural coordinates
    ap = f_align_lp_gc(ap,[apg(1,:),apg(end,:)],sctr);
    for seg = 1:length(ap)-1
      p = ap(seg:seg+1,:);
      pg = [apg(seg,:),apg(seg+1,:)];
      [W,Q] = quadrature(wall_int,'GAUSS',1) ;
      % find the distance between the two intersects (should be able to do this with det(J)
      [l,nv,mv,nnt,nmt,mmt] = f_segment_dist(pg);
      JO = l/2;
      mT = mel_elems(ii,3);
      kk = mel_elems(ii,1); 
      
      for k_in = 1:length(Q)
        [N1,dNdx1]=lagrange_basis('L2',Q(k_in));
        gpt = N1'*p;
        [N,dNdxi] = lagrange_basis(elemType,gpt) ;
        Nmat = enrNmat(N,iel,type_elem,enr_node(:,kk),elem_crk,xVertex,xTip,kk,modpen);
        gn = nv*Nmat*2*u(A);
        
        fn = Cm1(1,1)*(gn/mT);
        elem_force(seg,iel,2*k_in-1) = -1*fn;
        Gint(A) = Gint(A) + (fn)*W(k_in)*det(JO)*Nmat'*nv';
        Kglobal(A,A) = Kglobal(A,A) + 2*(Cm1(1,1)/mT)*W(k_in)*Nmat'*nnt*Nmat*det(JO) ;
        if elem_tan
          gt = mv*Nmat*2*u(A);
          ft = Cm1(1,2)*(gt/mT);
          elem_force(seg,iel,2*k_in) = -1*ft;
          Gint(A) = Gint(A) + (ft)*W(k_in)*det(JO)*Nmat'*mv';
          Kglobal(A,A) = Kglobal(A,A) + 2*(Cm1(1,2)/mT)*W(k_in)*Nmat'*mmt*Nmat*det(JO)/mT ;
        end
      end
      %disp(['iel is ',num2str(iel)]);
      %disp(['mT: ',num2str(mT)]);
      %disp(['gn: ',num2str(gn)]);
      %disp(['gt: ',num2str(gt)]);
    end
  end
end
%Plot Gauss points for checking
% plot mesh with crack and enriched nodes
% plotCrack(xCrk,enrich_node,plotmesh) ;
% plot(q(:,1),q(:,2),'r*') ;
% clear q
%out_str = ['The ratio of contact elements : ',num2str(cct/cc)];
%disp(out_str);
%fprintf(output_file,[out_str,'\n'])
