function [Kglobal,Fext] = KmatXFEM(enrich_node,elem_crk,type_elem,xTip,xVertex,...
    split_elem,tip_elem,vertex_elem,corner_elem,crack_node,enr_domain,pos,xCrk,Kglobal,Fext)

%declare global variables here
global node element numnode numelem elemType
global E C nu
global typeProblem typeCrack
global plotmesh plotNode
global gporder numtri
global plothelp
global orig_nn
global OPT rift_wall_pressure

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
[W2,Q2] = quadrature(2,intType,2) ;
% number of non-enriched df
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
            B = [B xfemBmat(Gpt,iel,type_elem,enrich_node(:,k),elem_crk,xVertex,crack_node,k)];
        end
        Ppoint =  N' * node(sctr,:);
        q = [q;Ppoint] ;
        try
          Kglobal(sctrB,sctrB) = Kglobal(sctrB,sctrB) + B'*C*B*W(kk)*det(JO) ;
        catch
          keyboard
        end

        switch OPT
        case 1
          if ~exist('ISSM_xx')
          sigma = f_extractStress(Ppoint)';
          end
          %keyboard
          Fext(sctrB) = Fext(sctrB) + B'*sigma*W(kk)*det(JO);
        case 2
          if ~exist('ISSM_xx')
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
              B2 = zeros(3,2*nn) ;
              B2(1,1:2:2*nn) = dNdx(:,1)' ;
              B2(2,2:2:2*nn) = dNdx(:,2)' ;
              B2(3,1:2:2*nn) = dNdx(:,2)' ;
              B2(3,2:2:2*nn) = dNdx(:,1)' ;
              Fext(sctrB(1:dfn)) = Fext(sctrB(1:dfn)) + B2'*sigma*W2(kk)*det(JO); % same jacobian
            end
          else
            Fext(sctrB) = Fext(sctrB) + B'*sigma*W(kk)*det(JO);
          end
        end
    end

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

%Plot Gauss points for checking
% plot mesh with crack and enriched nodes
% plotCrack(xCrk,enrich_node,plotmesh) ;
% plot(q(:,1),q(:,2),'r*') ;
% clear q
%if ~isempty(tan_element)
  %ntan = size(tan_element,1);
  %for iel = 1:ntan
      %sctr = tan_element(iel,:) ;
      %nn = length(sctr) ;
      %ke = 0 ;

      %%choose Gauss quadrature rules for elements
      %if strcmp(elemType,'Q4') 
        %intType = 'GAUSS' ;
      %elseif strcmp(elemType,'T3')
        %intType = 'TRIANGULAR' ;
      %end
      %[W,Q] = quadrature(IntOrder,intType,2) ;

      %%split_corner = ismember(iel,corner_elem) & ismember(iel,split_elem)

      %sctrB = [ ] ;
      %for k = 1:size(xCrk,2)
          %sctrB = [sctrB tan_assembly(iel,tan_element,pos(:,k))] ;
      %end

      %%loop over Gauss points
      %for kk = 1:size(W,1)
        %B = [] ;
        %Gpt = Q(kk,:) ;
        %[N,dNdxi] = lagrange_basis(elemType,Gpt) ;
        %JO = node(sctr,:)'*dNdxi ;
        %for k = 1:size(xCrk,2)
            %B = [B tan_xfemBmat(Gpt,iel,tan_elem,tan_elem_crk,crack_nodes)];
        %end
        %Kglobal(sctrB,sctrB) = Kglobal(sctrB,sctrB) + B'*C*B*W(kk)*det(JO) ;
        %if OPT == 1
          %sigma = f_extractStress(Ppoint)';
          %Fext(sctrB) = Fext(sctrB) + B'*sigma*W(kk)*det(JO);
        %end
      %end
  %end
%end
