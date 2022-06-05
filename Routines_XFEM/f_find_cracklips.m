function [ crack_lips, Flag_pen ] = f_find_cracklips( u, xCr, kk, enr_dom, type_elem, xCrl,xTip,xVertex,enr_node,crack_node,pos,split_elem, vertex_elem, tip_elem)
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Thu Nov 25 11:49:40 NZDT 2021
global node element elemType
global E nu C sigmato
global Jint iMethod
global epsilon

if strcmp(elemType,'Q4')
  corner = [1 2 3 4 1] ;
  nnode = [-1 -1;1 -1;1 1;-1 1] ;
elseif strcmp(elemType, 'T3')
  corner = [1 2 3 1] ;
  nnode = [0 0;1 0;0 1] ;
end

%Bfem = zeros(3,2*nn) ;
%Bfem(1,1:2:2*nn) = dNdx(:,1)' ;
%Bfem(2,2:2:2*nn) = dNdx(:,2)' ;
%Bfem(3,1:2:2*nn) = dNdx(:,2)' ;
%Bfem(3,2:2:2*nn) = dNdx(:,1)' ;
%BI_enr_p= [dNdx(in,1) 0 ; 0 dNdx(in,2) ; dNdx(in,2) dNdx(in,1)];
%BI_enr_n = [-1*dNdx(in,1) 0 ; 0 -1*dNdx(in,2) ; -1*dNdx(in,2) -1*dNdx(in,1)];

%num_cracks = size(xcr,2)
%num_enrich = size(enrdomain,1) 
elems = union(split_elem,vertex_elem);
elems = union(elems,tip_elem);

crack_lips = zeros( size(elems,1),6,4,size(xCr,2));
Flag_pen = 0;

for ii=1:length(elems)
  p = [];
  iel = elems(ii) ;
  sctr=element(iel,:);
  nn = length(sctr);
  vv = node(sctr,:);
  [phi] = dista(iel,xCrl) ;
  % First we find out what the points are that we are going to evalaute the crack lips at
  if ismember(iel, tip_elem)
    tip = xTip(iel,:);
    ntip = f_naturalpoint(tip,vv,20,1e-6);
    psi = f_dista2(iel,xCrl,tip);
    [cutEdge,nnodes] = f_edgedetect(nnode, corner,  phi, psi) ;
     p = [nnodes(end,:) , ntip ];
  else 
    [cutEdge, nnodes] = f_edgedetect(nnode, corner,  phi) ;
    nEdge = length(cutEdge);
    if nEdge==1 % then one of the nodes must be a crack_node
      crack_n = intersect(crack_node,sctr);
      if isempty(crack_n)
        error("only one edge but no crack node in elem ,",num2str(iel),", when evaluating crack lips")
      end
      crack_c = find(sctr==crack_n);
      nnodes = [ nnodes; nnodes(crack_c,:) ] ;
    end

    p = [nnodes(end-1,:);nnodes(end,:)];

    if ismember(iel,vertex_elem)
      tip = xVertex(iel,:);
      ntip = f_naturalpoint(tip,vv,20,1e-6);
      p = [p(1,:) ; ntip ; p(2,:) ];
    end
  end

  % now we are going to evaluate the displacement of each point for the positive, negative and midpoint of the crack
 [A,BrI,QT,Tip,alpha] = f_enrich_assembly(iel,pos,type_elem,elem_crk,enr_node);


  for gp = 1:size(p,1)
    c_inds = [2*gp - 1, 2*gp ];
    gpt = p(gp,:);
    [N,dNdxi] = lagrange_basis(elemType,gpt);
    cc_m = N'*node(sctr,:);
    if any(enr_node(sctr)==1)
      up_down = 0;
      in = find(enr_node(sctr)==1,1);
      if type_elem(iel,1) == 1   %looking for the "tip" element
        ref_elem = iel;
        Rpt = 1;
      else    %trovo l'elemento/fessura a cui fa riferimento il nodo (SOLO 1 RIF AUTORIZZATO!!)
        [sctrn,xx] = find(element == sctr(in));
        [ele,xx] = find(type_elem(sctrn,:)==1);
        ref_elem = sctrn(ele);
        elem_blend = 1;
        nR = find(enr_node(sctr)==1);
        Rpt = sum(N(nR));
      end
      % compute branch functions at Gauss point
      xCre  = [xCrl(ref_elem,1) xCrl(ref_elem,2); xCrl(ref_elem,3) xCrl(ref_elem,4)];
      seg   = xCre(2,:) - xCre(1,:);
      alpha = atan2(seg(2),seg(1));
      Tip  = [xCre(2,1) xCre(2,2)];
      QT    = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
      xp    = QT*(cc_m-Tip)';           % crack aligned coordinates
      r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2)); %could almost use xp(2) here but might not work when kinked close to tip 
      theta = atan2(xp(2),xp(1));

      
      if ( theta > pi | theta < -pi)
          disp (['something wrong with angle ',num2str(thet)]);
      end
      if abs(abs(theta) - pi) < 0.001
       [Br_u,dBdx,dBdy] = branch_gp(r,pi,alpha);
       [Br_d,dBdx,dBdy] = branch_gp(r,-1*pi,alpha);
       up_down = 1;
      else 
       [Br_u,dBdx,dBdy] = branch_gp(r,theta,alpha);
       Br_d = Br_u;
      end
    end

    [N,dNdxi] = lagrange_basis(elemType,gpt);
    Ue = element_disp(iel,pos(:,kk),enr_node(:,kk),u,kk);
    c_m = [N'*Ue(1:2:2*nn),N'*Ue(2:2:2*nn)]; %u_fem
    c_p = [ 0, 0];
    c_d = [ 0, 0];
    idx = nn+1;
    for in = 1:nn
      nodeI = sctr(in) ;
      Ni = N(in);
      if (enr_node(nodeI) == 2) | (enr_node(nodeI) == 3)    % H(x) enriched node
          c_p = c_p + Ni*[Ue(idx*2-1),Ue(idx*2)];
          c_d = c_d + 0*Ni*[Ue(idx*2-1),Ue(idx*2)];
          idx = idx + 1;
      elseif(enr_node(nodeI) == 1)
        xp    = QT*(node(sctr(in),:)-Tip)';
        r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
        theta = atan2(xp(2),xp(1));
        
        if ( theta > pi | theta < -pi)
            disp (['something wrong with angle ',num2str(thet)]);
        end
        [BrI] = branch_node(r,theta);

        for i = 1:4
          c_p = c_p + Rpt*Ni*(Br_u(i)-BrI(i)) * [Ue(idx*2-1),Ue(idx*2)];
          c_d = c_d + Rpt*Ni*(Br_d(i)-BrI(i)) * [Ue(idx*2-1),Ue(idx*2)];
          idx = idx + 1;
          %if (abs(cc_m(1)-1.) < 0.001)
            %disp(['element ',num2str(iel),'   function ',num2str(i),' :'])
            %disp(['Br = ', num2str(Br(i))])
          %end  
        end
      end
    end
      c_plus = c_m + c_p;
      c_down = c_m + c_d; 
      crack_lips(ii,c_inds,1,kk) = gpt ;
      crack_lips(ii,c_inds,2,kk) = c_plus;
      crack_lips(ii,c_inds,3,kk) = c_down;
      crack_lips(ii,c_inds,4,kk) = cc_m;
      % check that there is no interpenetration
      [phi] = f_dista_point( cc_m+c_p, iel, xCrl, 1e-8 ); 
      gn = c_p - c_d;
      if any(gn < 0) % crack tips can be a problem for this 
        Flag_pen = 1;
      end
    end
end
