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
elems = union(elems,tip_elem);
if melangeforce
  elemst = tan_elem;
  mel_elems = []

  for kk = 1:size(xCrk,2)

    for iel=1:length(elems)                     %loop on elems (=elements selected for enrichment)
      found = 0;
      for kj = 1:size(xM.coor,1)-1       %loop over the elements of the fracture
        e = elems(iel);
        if xM.melange(kj)
          q1 = xM(kk).coor(kj,:); 
          q2 = xM(kk).coor(kj+1,:);
          [flag1,flag2,cn_] = crack_interact_element([q1,q2],e,[]);
          if flag1
            if found 
              mel_elems(end,2) = (mel_elems(end,2) + xM.width(kj))/2
              break
            else
              found = 1;
              mel_elems = [mel_elems; e xM.width(kj) ];
            end
          elseif found % the only chance of an element belonging two 2 crack sections is if they are one after the other 
             break
          end
        end
      end
    end
    %for iel=1:length(elemst)                     %loop on elems (=elements selected for enrichment)
      %for kj = 1:size(xM.coor,1)-1       %loop over the elements of the fracture
        %if xM.melange(kj)
          %e = elemst(iel);
          %q1 = xM(kk).coor(kj,:); 
          %q2 = xM(kk).coor(kj+1,:);
          %[flag1,flag2,cn_] = crack_interact_element([q1,q2],e,[]);
          %if flag2
              %mel_elems = [mel_elems; e xM.width(kj)/2 ];
              %break
          %end
        %end
      %end
    %end
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
      cc = cc + 1;
      iel = elems(ii) ;
      sctr = element(iel,:) ;
      skip = 0;
      nn = length(sctr) ;

      [A,BrI,QT,Tip,alpha] = f_enrich_assembly(iel,pos,type_elem,elem_crk,enr_node);
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
          Nmat = enrNmat(N,iel,type_elem,enr_node(:,kk),elem_crk,xVertex,kk,true);
          gn = nv*Nmat*2*u(A);
          if gn < 0
            if k_in ==1
              if seg == 1
                cct = cct + 1;
              end
            end
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
end
if melangeforce
  for kk = 1:size(xCrk,2) %what's the crack?
    for ii=1:size(mel_elems,1)
      cc = cc + 1;
      iel = mel_elems(ii) ;
      sctr = element(iel,:) ;
      skip = 0;
      nn = length(sctr) ;

      [A,BrI,QT,Tip,alpha] = f_enrich_assembly(iel,pos,type_elem,elem_crk,enr_node);
      [ap,apg] = f_crack_wall(iel,nnode,corner,tip_elem,vertex_elem,elem_crk,xTip,xVertex,crack_node); % elem_crk in natural coordinates
      for seg = 1:length(ap)-1
        p = ap(seg:seg+1,:);
        pg = [apg(seg,:),apg(seg+1,:)];
        [W,Q] = quadrature(wall_int,'GAUSS',1) ;
        % find the distance between the two intersects (should be able to do this with det(J)
        [l,nv,mv,nnt,nmt,mmt] = f_segment_dist(pg);
        JO = l/2;


        if melangeforce
          if ismember(iel,mel_elems(:,1))
            ind = find(mel_elems(:,1)==iel);
            mT = mel_elems(ind,2);
          else
            continue
          end
          
          for k_in = 1:length(Q)
            gpt = gpts(k_in,:) ;
            [N,dNdxi] = lagrange_basis(elemType,gpt) ;
            Nmat = [N(1), 0, N(2), 0, N(3), 0 ; 0, N(1), 0 , N(2), 0, N(3)];
            Nmat = enrNmat(N,iel,type_elem,enr_node(:,kk),elem_crk,xVertex,kk);
            gn = nv*Nmat*2*u(A);
            gt = mv*Nmat*2*u(A);
            fn = Cm1(1,1)*(gn/mT);
            ft = Cm1(1,2)*(gt/mT);
            elem_force(seg,iel,2*k_in-1) = Cm1(1,1)*gn;
            elem_force(seg,iel,2*k_in) = Cm1(1,2)*gt;
            Kglobal(A,A) = Kglobal(A,A) + Cm1(1,1)*W(k_in)*Nmat'*nnt*Nmat*det(JO)/mT ;
            Kglobal(A,A) = Kglobal(A,A) + Cm1(1,2)*W(k_in)*Nmat'*mmt*Nmat*det(JO)/mT ;
            Gint(A) = Gint(A) + W(k_in)*det(JO)*fn*Nmat'*nv';
            Gint(A) = Gint(A) + W(k_in)*det(JO)*ft*Nmat'*mv';
          end
          disp(['iel is ',num2str(iel)]);
          disp(['mT: ',num2str(mT)]);
          disp(['gn: ',num2str(gn)]);
          disp(['gt: ',num2str(gt)]);
        end
      end
    end
  end
end
%Plot Gauss points for checking
% plot mesh with crack and enriched nodes
% plotCrack(xCrk,enrich_node,plotmesh) ;
% plot(q(:,1),q(:,2),'r*') ;
% clear q
ratio_c = cct/cc;
%out_str = ['The ratio of contact elements : ',num2str(cct/cc)];
%disp(out_str);
%fprintf(output_file,[out_str,'\n'])
