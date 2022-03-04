function [Kglobal,Gint] = KTmatXFEM(E_pen,enr_node,elem_crk,type_elem,xTip,xVertex,...
    split_elem,tip_elem,vertex_elem,corner_elem,crack_nodes,pos,xCrk,Kglobal,u)

%declare global variables here
global node element numnode numelem elemType
global E C nu
global typeProblem typeCrack
global plotmesh plotNode
global gporder numtri
global plothelp
global orig_nn

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
%loop over enriched elements
elems = union(split_elem,vertex_elem);
% Gint
Gint =  zeros(size(u))

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
    sctrA = [];
    skip = 0;
    nn = length(sctr) ;
    for ni = 1:nn
      nodeI = sctr(ni); 
      if (enr_node(nodeI)==1) | (enr_node(nodeI)==0)
        skip = 1;
      end
      AA = [2*pos(nodeI)-1;2*pos(nodeI)];
      sctrA = [sctrA; AA];
    end
    if skip
      if plothelp
         disp(['Element Skipped for penalty: different enrichments found iel = ',num2str(iel)])
      end
      continue
    end


    vv = node(sctr,:);
    [phi] = dista(iel,elem_crk) ;
    if ismember(iel, tip_elem) % for now we wont deal with this element
      tip = xTip(iel,:);
      ntip = f_naturalpoint(tip,vv,20,epsilon);
      psi = f_dista2(iel,elem_crk,tip);
      [cutEdge,nnodes] = f_edgedetect(nnode, corner,  phi, psi) ;
      p = [nnodes(end,:); ntip ];
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

      if ismember(iel,vertex_elem) % we are jsut going to draw a ligne across the two intersections
        tip = xVertex(iel,:);
        ntip = f_naturalpoint(tip,vv,20,epsilon);
        %p = [p(1,:) ; ntip ; p(2,:) ];
      end
    end

    % need to find the normal
    

    [W,Q] = quadrature(2,'GAUSS',1) ;
    [N1,dNdx1]=lagrange_basis('L2',Q(1));
    [N2,dNdx2]=lagrange_basis('L2',Q(2));
    gpts = [N1'*p; N2'*p];
    % find the distance between the two intersects (should be able to do this with det(J)
    x0 = elem_crk(iel,1) ; y0 = elem_crk(iel,2) ;
    x1 = elem_crk(iel,3) ; y1 = elem_crk(iel,4) ;
    l = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0)) ;
    nv = [(y0-y1),(x1-x0)]./l;
    nnt = nv'*nv;

    JO = l/2; % check this isn't in natural

    for kk = 1:2
        gpt = gpts(kk,:) ;
        [N,dNdxi] = lagrange_basis(elemType,gpt) ;
        Nmat = [N(1), 0, N(2), 0, N(3), 0 ; 0, N(1), 0 , N(2), 0, N(3)];
        gn = nv*Nmat*u(sctrA)
        Ppoint =  N' * node(sctr,:);
        if gn < 0
        try
          Kglobal(sctrA,sctrA) = Kglobal(sctrA,sctrA) + E_pen*W(kk)*Nmat'*nnt*Nmat*det(JO) ;
          Gint(sctrA) = Gint(sctrA) + E_pen*W(kk)*det(JO)*gn*Nmat'*nv';
        catch
          keyboard
        end
        end

        
      if plothelp
      figure(2)
      %plotMesh(node,element(iel,:),'T3','r-','no')
      %if ismember(iel,split_elem)   

        %ppl = plot(Ppoint(1),Ppoint(2),'*m','linestyle','none','markersize',2)
        %%keyboard
        %delete(ppl)
      %end
      plot(Ppoint(1),Ppoint(2),'*c','linestyle','none','markersize',1)
      keyboard
      
    end
    end
  end
end
%Plot Gauss points for checking
% plot mesh with crack and enriched nodes
% plotCrack(xCrk,enrich_node,plotmesh) ;
% plot(q(:,1),q(:,2),'r*') ;
% clear q
