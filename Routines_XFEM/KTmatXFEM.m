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

for kk = 1:size(xCrk,2) %what's the crack?
  for ii=1:size(elems,1)
    %iel = elems(ii) ;
    %sctr=element(iel,:);
    %nn = length(sctr);
    %p1 = crack_lips(ii,1:2,4,kk);
    %p2 = crack_lips(ii,3:4,4,kk);
    %p3 = crack_lips(ii,5:6,4,kk);
    %n1 = crack_lips(ii,1:2,1,kk);
    %n2 = crack_lips(ii,3:4,1,kk);
    %n3 = crack_lips(ii,5:6,1,kk);
    %dup1 = crack_lips(ii,1:2,2,kk);
    %dup2 = crack_lips(ii,3:4,2,kk);
    %dup3 = crack_lips(ii,5:6,2,kk);
    %ddown1 = crack_lips(ii,1:2,3,kk);
    %ddown2 = crack_lips(ii,3:4,3,kk);
    %ddown3 = crack_lips(ii,5:6,3,kk);

    %[fd_xy,fu_xy] = f_crack_contact_penalty([p1,p2],[ddown1,ddown2], [dup1,dup2],F_app);
    %Fc = f_apply_crack_force(Fc,fd_xy,fu_xy,n1,n2,iel,xCr,xCrl,xTip,pos,type_elem,enr_node);


    %if ismember(iel,vertex_elem)
      %[fd_xy,fu_xy] = f_crack_contact_penalty([p2,p3],[ddown2,ddown3], [dup2, dup3], F_app);
      %Fc = f_apply_crack_force(Fc,fd_xy,fu_xy,n2,n3,iel,xCr,xCrl,xTip,pos,type_elem,enr_node);

    %end
    iel = elems(ii) ;
    sctr = element(iel,:) ;
    skip = 0;
    nn = length(sctr) ;

    [A,BrI,QT,Tip,alpha] = f_enrich_assembly(iel,pos,type_elem,elem_crk,enr_node);
    p = f_crack_wall(iel,nnode,corner,tip_elem,vertex_elem,elem_crk,xTip,crack_node); % elem_crk in natural coordinates
    %vv = node(sctr,:);
    %[phi] = dista(iel,elem_crk) ;
    %if ismember(iel, tip_elem) % for now we wont deal with this element
      %tip = xTip(iel,:);
      %ntip = f_naturalpoint(tip,vv,20,epsilon);
      %psi = f_dista2(iel,elem_crk,tip);
      %[cutEdge,nnodes] = f_edgedetect(nnode, corner,  phi, psi) ;
      %p = [nnodes(end,:); ntip ];
    %else 
      %[cutEdge, nnodes] = f_edgedetect(nnode, corner,  phi) ;
      %nEdge = length(cutEdge);
      %if nEdge==1 % then one of the nodes must be a crack_node
        %crack_n = intersect(crack_node,sctr);
        %if isempty(crack_n)
          %error("only one edge but no crack node in elem ,",num2str(iel),", when evaluating crack lips")
        %end
        %crack_c = find(sctr==crack_n);
        %nnodes = [ nnodes; nnodes(crack_c,:) ] ;
      %end

      %p = [nnodes(end-1,:);nnodes(end,:)];

      %if ismember(iel,vertex_elem) % we are jsut going to draw a ligne across the two intersections
        %tip = xVertex(iel,:);
        %ntip = f_naturalpoint(tip,vv,20,epsilon);
        %%p = [p(1,:) ; ntip ; p(2,:) ];
      %end
    %end

    % need to find the normal
     
    for seg = 1:length(p)-1

      [W,Q] = quadrature(3,'GAUSS',1) ;
      % find the distance between the two intersects (should be able to do this with det(J)
      [l,nv,mv,nnt,nmt,mmt] = f_segment_dist(elem_crk(iel,:));
      JO = l/2;

      if contact 
        for k_in = 1:length(Q)
          [N1,dNdx1]=lagrange_basis('L2',Q(k_in));
          gpt = N1'*p;
      % find the distance between the two intersects (should be able to do this with det(J)
          [N,dNdxi] = lagrange_basis(elemType,gpt) ;
          pint =  N' * node(sctr,:);
          Nmat = enrNmat(N,iel,type_elem,enr_node(:,kk),elem_crk,xVertex,kk,true)
          gn = nv*Nmat*2*u(A);
          if gn < 0
            elem_force(seg,iel,2*k_in-1) = -1*E_pen*gn;
            Gint(A) = Gint(A) + E_pen*W(k_in)*det(JO)*gn*Nmat'*nv';
            Kglobal(A,A) = Kglobal(A,A) + 2*E_pen*W(k_in)*Nmat'*nnt*Nmat*det(JO) ;
            % stabalization term
            %Kglobal(A,A) = Kglobal(A,A) - W(k_in)*((E_pen^2)/(2*E))*((Nmat'-1/3)*nnt*(Nmat-1/3))*det(JO);
            %if any(enrich_node(sctr)==1)
              %xp    = QT*(pint-Tip)';           % local coordinates
              %r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
              %theta = atan2(xp(2),xp(1));
              
              %if ( theta > pi | theta < -pi)
                  %disp (['something wrong with angle ',num2str(thet)]);
              %end
              %if abs(abs(theta) - pi) < 0.001
               %[Br_u,dBdx,dbdy] = branch_gp(r,pi,alpha);
               %[Br_d,dBdx,dbdy] = branch_gp(r,-1*pi,alpha);
              %else 
               %[Br_u,dbdx,dbdy] = branch_gp(r,theta,alpha);
               %Br_d = Br_u;
              %end
            %end
            %nA = [1,2];
            %for ni = 1:nn
              %if n1(ni)
                %n_row = sum(n1(1:ni));
                %for i = 1:4
                  %Kglobal(A(nA),A(nA)) = Kglobal(A(nA),A(nA)) + E_pen*W(k_in)*N(ni)*nnt*N(ni)*det(JO) ;
                  %Gint(A(nA)) = Gint(A(nA)) + E_pen*W(k_in)*det(JO)*gn*N(ni)*nv';
                  %nA = [nA(1)+2,nA(2)+2];
                %end
                %elem_force(iel,2*k_in-1) = E_pen*gn;
              %else
                %elem_force(iel,2*k_in-1) = E_pen*gn;
                %Gint(A(nA)) = Gint(A(nA)) + E_pen*W(k_in)*det(JO)*gn*N(ni)*nv';
                %nA = [nA(1)+2,nA(2)+2];
              %end 
            %end
          end
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

      if melangeforce
        if ismember(iel,mel_elems(:,1))
          ind = find(mel_elems(:,1)==iel);
          mT = mel_elems(ind,2);
        else
          continue
        end
        keyboard
        
        for k_in = 1:2
          gpt = gpts(k_in,:) ;
          [N,dNdxi] = lagrange_basis(elemType,gpt) ;
          Nmat = [N(1), 0, N(2), 0, N(3), 0 ; 0, N(1), 0 , N(2), 0, N(3)];
          Nmat = enrNmat(N,iel,type_elem,enr_node(:,kk),elem_crk,xVertex,kk)
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
%Plot Gauss points for checking
% plot mesh with crack and enriched nodes
% plotCrack(xCrk,enrich_node,plotmesh) ;
% plot(q(:,1),q(:,2),'r*') ;
% clear q
