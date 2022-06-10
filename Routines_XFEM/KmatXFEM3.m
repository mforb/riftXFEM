function [Kglobal,Fext] = KmatXFEM3(enrich_node,elem_crk,type_elem,xTip,xVertex,...
    split_elem,tip_elem,vertex_elem,corner_elem,crack_node,enr_domain,pos,xCrk,Kglobal,Fext)

%declare global variables here
global node element numnode numelem elemType
global E C nu
global typeProblem typeCrack
global plotmesh plotNode
global gporder numtri
global plothelp
global orig_nn
global OPT rift_wall_pressure ISSM_xx

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
% for mapping options 2 and 3 we are going to create some variables to treat elements as if they are not enriched
%enrich_node2 = zeros(size(enrichNode));
    % for import option 2
if strcmp(elemType,'Q4') 
  intType = 'GAUSS' ;
  corner = [1 2 3 4 1] ;
  nnode = [-1 -1;1 -1;1 1;-1 1] ;
else
  intType = 'TRIANGULAR';
  corner = [1 2 3 1] ;
  nnode = [0 0;1 0;0 1] ;
end
% number of non-enriched df
[W2,Q2] = quadrature(2,intType,2) ;
dfn = size(W2,1)*2;


q = [] ;
%loop over ALL elements
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
        sctrB = [sctrB assembly(iel,enrich_node(:,k),pos(:,k),k,crack_node)] ;
    end

    %loop over Gauss points
    for kk = 1:size(W,1)
        B = [] ;
        Gpt = Q(kk,:) ;
        [N,dNdxi] = lagrange_basis(elemType,Gpt) ;
        JO = node(sctr,:)'*dNdxi ;
        for k = 1:size(xCrk,2)
            B = [B xfemBmat(Gpt,iel,type_elem,enrich_node(:,k),elem_crk,xVertex,xTip,crack_node,k)];
        end
        Ppoint =  N' * node(sctr,:);
        q = [q;Ppoint] ;
        Kglobal(sctrB,sctrB) = Kglobal(sctrB,sctrB) + B'*C*B*W(kk)*det(JO) ;

        switch OPT
        case 1
          if ~isempty('ISSM_xx')
            sigma = f_getstress(iel);
          else
            sigma = f_extractStress(Ppoint)';
          end
          %keyboard
          Fext(sctrB) = Fext(sctrB) + B'*sigma*W(kk)*det(JO);
        case 2
          if ~isempty('ISSM_xx')
            sigma = f_getstress(iel);
          else
            sigma = f_extractStress(Ppoint)';
          end
          if ismember(iel,enr_domain)  % over kill, maybe should check if any nodes are enriched
            % pretend the enrdomain is not enriched
            if kk <= dfn/2
              Gpt = Q2(kk,:);
              [N,dNdxi] = lagrange_basis(elemType,Gpt);  % element shape functions
              %J0 = node(sctr,:)'*dNdxi;                 % element Jacobian matrix
              invJO = inv(JO);
              dNdx  = dNdxi*invJO;                      % derivatives of N w.r.t XY
              B2 = B(:,1:2*nn);
              sctrB2 = sctrB(1:2*nn);
              Fext(sctrB2) = Fext(sctrB2) + B2'*sigma*W2(kk)*det(JO); % same jacobian
            end
          else
            Fext(sctrB) = Fext(sctrB) + B'*sigma*W(kk)*det(JO);
          end
        end
    end

end

