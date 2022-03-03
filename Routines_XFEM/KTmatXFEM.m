function [KTglobal] = KTmatXFEM(K,enrich_node,elem_crk,type_elem,xTip,xVertex,...
    split_elem,tip_elem,vertex_elem,corner_elem,crack_nodes,pos,xCrk,Kglobal)

%declare global variables here
global node element numnode numelem elemType
global E C nu
global typeProblem typeCrack
global plotmesh plotNode
global gporder numtri
global plothelp
global orig_nn

% we are adding the gradient dt/du to K
KTglobal = K

%loop over enriched elements
elems = union(split_elem,vertex_elem);

for kk = 1:size(xCr,2) %what's the crack?
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
    sctrB = []
    for nodeI = 1:length(sctr) 
      AA = [u(2*pos(nodeI)-1);u(2*pos(nodeI))];
      sctrA = [sctrA; AA]
    end
    nn = length(sctr) ;
    vv = node(sctr,:);
    [phi] = dista(iel,xCrl) ;
    if ismember(iel, tip_elem) % for now we wont deal with this element
      tip = xTip(iel,:);
      ntip = f_naturalpoint(tip,vv,20,epsilon);
      psi = f_dista2(iel,xCrl,tip);
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
    gpts = [N1'*p; N2'*p]
    % find the distance between the two intersects (should be able to do this with det(J)
    x0 = elem_crk(iel,1) ; y0 = elem_crk(iel,2) ;
    x1 = elem_crk(iel,3) ; y1 = elem_crk(iel,4) ;
    nv = [(y0-y1),(x1-x0)]./l
    nnt = n'*n

    l = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0)) ;
    JO = l/2; % check this isn't in natural

    for kk = 1:2
        gpt = gpts(kk,:) ;
        [N,dNdxi] = lagrange_basis(elemType,gpt) ;
        Nmat = [N(1), 0, N(2), 0, N(3), 0 ; 0, N(1), 0 , N(2), 0, N(3)];
        try
          Kglobal(sctrA,sctrA) = Kglobal(sctrA,sctrA) + W(kk)*Nmat'*nnt*Nmat*det(JO) ;
        catch
          keyboard
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
      
    end
    end
  end
end
%Plot Gauss points for checking
% plot mesh with crack and enriched nodes
% plotCrack(xCrk,enrich_node,plotmesh) ;
% plot(q(:,1),q(:,2),'r*') ;
% clear q
