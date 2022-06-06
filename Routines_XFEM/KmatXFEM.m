function [Kglobal] = KmatXFEM(enrich_node,elem_crk,type_elem,xTip,xVertex,...
    split_elem,tip_elem,vertex_elem,corner_elem,crack_nodes,pos,xCrk,Kglobal)

%declare global variables here
global node element numnode numelem elemType
global E C nu
global typeProblem typeCrack
global plotmesh plotNode
global gporder numtri
global plothelp
global orig_nn

if plothelp
  figure(2)
  hold on
  split_nodes = find(enrich_node == 2);
  tip_nodes   = find(enrich_node == 1);
  n1 = plot(node(split_nodes,1),node(split_nodes,2),'r*');
  n2 = plot(node(tip_nodes,1),node(tip_nodes,2),'rs');
  set(n1,'MarkerSize',5);
  set(n2,'MarkerSize',5);
  plotMesh_numbered(node,element,elemType,'b-','no')
end

q = [] ;
%loop over elements
for iel = 1:numelem
    sctr = element(iel,:) ;
    nn = length(sctr) ;
    ke = 0 ;

    %choose Gauss quadrature rules for elements
    [W,Q] = gauss_rule(iel,enrich_node,elem_crk,...
        xTip,xVertex,tip_elem,split_elem,vertex_elem,corner_elem,xCrk) ;

    %split_corner = ismember(iel,corner_elem) & ismember(iel,split_elem)

    sctrB = [ ] ;
    for k = 1:size(xCrk,2)
        sctrB = [sctrB assembly(iel,enrich_node(:,k),pos(:,k),k,crack_nodes)] ;
    end

    %loop over Gauss points
    for kk = 1:size(W,1)
        B = [] ;
        Gpt = Q(kk,:) ;
        [N,dNdxi] = lagrange_basis(elemType,Gpt) ;
        JO = node(sctr,:)'*dNdxi ;
        for k = 1:size(xCrk,2)
            B = [B xfemBmat(Gpt,iel,type_elem,enrich_node(:,k),elem_crk,xVertex,xTip,crack_nodes,k)];
        end
        Ppoint =  N' * node(sctr,:);
        q = [q;Ppoint] ;
        Kglobal(sctrB,sctrB) = Kglobal(sctrB,sctrB) + B'*C*B*W(kk)*det(JO) ;

        %if ismember(iel,[106,246,184])
          %nn   = length(sctr);
          %cnt = 0 ;

          %for k = 1 : nn
            %cnt = cnt + 1 ;
            %sctrB_u(2*k-1) = 2*sctr(k)-1 ;
            %sctrB_u(2*k)   = 2*sctr(k)   ;
            %sctrB_a(2*cnt - 1) = 2 * pos(sctr(k)) - 1;
            %sctrB_a(2*cnt    ) = 2 * pos(sctr(k))    ;
          %end
          %[B_u,B_a] = tan_xfemBmat(Gpt,iel,element,elem_crk,crack_nodes)
          %Kglobal(sctrB_u,sctrB_u) = Kglobal(sctrB_u,sctrB_u) + B_u'*C*B_u*W(kk)*det(JO) ;
          %Kglobal(sctrB_a,sctrB_u) = Kglobal(sctrB_a,sctrB_u) + B_a'*C*B_u*W(kk)*det(JO) ;
          %Kglobal(sctrB_u,sctrB_a) = Kglobal(sctrB_u,sctrB_a) + B_u'*C*B_a*W(kk)*det(JO) ;
          %Kglobal(sctrB_a,sctrB_a) = Kglobal(sctrB_a,sctrB_a) + B_a'*C*B_a*W(kk)*det(JO) ;
          %Kglobal2(sctrB,sctrB) = Kglobal2(sctrB,sctrB) + B'*C*B*W(kk)*det(JO) ;
        %else
           %Kglobal(sctrB,sctrB) = Kglobal(sctrB,sctrB) + B'*C*B*W(kk)*det(JO) ;
           %Kglobal2(sctrB,sctrB) = Kglobal2(sctrB,sctrB) + B'*C*B*W(kk)*det(JO) ;
        %end

        
      %if plothelp
      %figure(2)
      %%plotMesh(node,element(iel,:),'T3','r-','no')
      %%if ismember(iel,split_elem)   

        %%ppl = plot(Ppoint(1),Ppoint(2),'*m','linestyle','none','markersize',2)
        %%%keyboard
        %%delete(ppl)
      %%end
      %plot(Ppoint(1),Ppoint(2),'*c','linestyle','none','markersize',1)
      
    %end
    end
end


% plot mesh with crack and enriched nodes
% plotCrack(xCrk,enrich_node,plotmesh) ;
% plot(q(:,1),q(:,2),'r*') ;
% clear q
