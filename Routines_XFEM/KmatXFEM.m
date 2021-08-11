function [Kglobal] = KmatXFEM(enrich_node,elem_crk,type_elem,xTip,xVertex,...
    split_elem,tip_elem,vertex_elem,pos,xCrk,Kglobal)

%declare global variables here
global node element numnode numelem elemType
global E C nu
global typeProblem typeCrack
global plotmesh plotNode
global gporder numtri

q = [] ;
%loop over elements
for iel = 1:numelem
    sctr = element(iel,:) ;
    nn = length(sctr) ;
    ke = 0 ;

    %choose Gauss quadrature rules for elements
    [W,Q] = gauss_rule(iel,enrich_node,elem_crk,...
        xTip,xVertex,tip_elem,split_elem,vertex_elem,xCrk) ;

    %Transfrom these Gauss points to global coords for plotting ONLY
    for igp = 1:size(W,1)
        gpnt = Q(igp,:) ;
        [N,dNdxi] = lagrange_basis(elemType,gpnt) ;
        Gpnt = N'*node(sctr,:) ;
        q = [q;Gpnt] ;
    end

    sctrB = [ ] ;
    for k = 1:size(xCrk,2)
        sctrB = [sctrB assembly(iel,enrich_node(:,k),pos(:,k),k)] ;
    end

    %loop over Gauss points
    for kk = 1:size(W,1)
        B = [] ;
        Gpt = Q(kk,:) ;
        [N,dNdxi] = lagrange_basis(elemType,Gpt) ;
        JO = node(sctr,:)'*dNdxi ;
        for k = 1:size(xCrk,2)
            B = [B xfemBmat(Gpt,iel,type_elem,enrich_node(:,k),elem_crk,xVertex,k)] ;
        end
        Ppoint =  N' * node(sctr,:);
            
        Kglobal(sctrB,sctrB) = Kglobal(sctrB,sctrB) + B'*C*B*W(kk)*det(JO) ;
    end
end

%Plot Gauss points for checking
% plot mesh with crack and enriched nodes
% plotCrack(xCrk,enrich_node,plotmesh) ;
% plot(q(:,1),q(:,2),'r*') ;
% clear q