function [f] = ForceVector(f)

%--- Purpose
%- compute nodal force vector

global sigmato E nu C P
global elemType
global typeProblem typeCrack
global node element numnode numelem bcNodes edgNodes

if( strcmp(typeProblem,'eCrkTen') || strcmp(typeProblem,'cCrkTen') ...
        || strcmp(typeProblem,'inCrkTen') || strcmp(typeProblem,'dEdCrkTen') ...
        || strcmp(typeProblem,'uSrDefined') || strcmp(typeProblem,'curvedCrk') )
    topEdge = edgNodes{3} ;
    [W,Q] = quadrature(1,'GAUSS',1) ;
    for e = 1:size(topEdge,1)
        sctr = topEdge(e,:) ;
        sctry = sctr.*2 ;
        for q = 1:size(W,1)
            pt = Q(q,:) ; wt = W(q) ;
            N = lagrange_basis('L2',pt) ;
            J0 = abs( node(sctr(2)) - node(sctr(1)) )/2 ;
            f(sctry) = f(sctry) + N*sigmato*det(J0)*wt ;
        end
    end
elseif( strcmp(typeProblem,'eCrkTen2') )
    topEdge = edgNodes{3} ;
    [W,Q] = quadrature(1,'GAUSS',1) ;
    for e = 1:size(topEdge,1)
        sctr = topEdge(e,:) ;
        sctry = sctr.*2 ;
        for q = 1:size(W,1)
            pt = Q(q,:) ; wt = W(q) ;
            N = lagrange_basis('L2',pt) ;
            J0 = abs( node(sctr(2)) - node(sctr(1)) )/2 ;
            f(sctry) = f(sctry) + N*sigmato*det(J0)*wt ;
        end
    end
elseif( strcmp(typeProblem,'Griffith') || strcmp(typeProblem,'dispfriction'))
    f = f ;
elseif( strcmp(typeProblem,'centre')  || strcmp(typeProblem,'centreF') )
    botEdge = edgNodes{1} ;
    [W,Q] = quadrature(1,'GAUSS',1) ;
    for e = 1:size(botEdge,1)
        sctr = botEdge(e,:) ;
        sctry = sctr.*2;
        for q = 1:size(W,1)
            pt = Q(q,:) ; wt = W(q) ;
            N = lagrange_basis('L2',pt) ;
            J0 = abs( node(sctr(2)) - node(sctr(1)) )/2 ;
            f(sctry) = f(sctry) - N*sigmato*det(J0)*wt ;
        end
    end
    if( strcmp(typeProblem,'centreF') )
      [W,Q] = quadrature(1,'TRIANGULAR',1) ;
      els = [600,5019,5044,5259,5459]
      for i = 1:length(els)
          e = els(i) ;
          sctr = element(e,:) ;
          sctry = sctr.*2;
          for q = 1:size(W,1)
              pt = Q(q,:) ; wt = W(q) ;
              [N,dNdxi] = lagrange_basis(elemType,pt);  % element shape functions
              J0 = node(sctr,:)'*dNdxi;                 % element Jacobian matrix
              f(sctry) = f(sctry) - .1*N*sigmato*det(J0)*wt ;
          end
      end
    end
elseif( strcmp(typeProblem,'eCrkShear') )
    topEdge = edgNodes{3} ;
    [W,Q] = quadrature(1,'GAUSS',1) ;
    for e = 1:size(topEdge,1)
        sctr = topEdge(e,:) ;
        sctry = sctr.*2 - 1;
        for q = 1:size(W,1)
            pt = Q(q,:) ; wt = W(q) ;
            N = lagrange_basis('L2',pt) ;
            J0 = abs( node(sctr(2)) - node(sctr(1)) )/2 ;
            f(sctry) = f(sctry) + N*sigmato*det(J0)*wt ;
        end
    end
elseif( strcmp(typeProblem,'ISSM') )
    %frontEdge = edgNodes{1} ;
    %[W,Q] = quadrature(1,'GAUSS',1) ;
    %for e = 1:size(frontEdge,1)
        %iel = f_findElemfromEdge(frontEdge(e,:))
        %sctr = frontEdge(e,:) ;
        %sctry = sctr.*2 - 1;
        %fh = f_getHeightF(iel);
        %for q = 1:size(W,1)
            %pt = Q(q,:) ; wt = W(q) ;
            %N = lagrange_basis('L2',pt) ;
            %J0 = abs( node(sctr(2)) - node(sctr(1)) )/2 ;
            %f(sctry) = f(sctry) + N*fh*det(J0)*wt ;
        %end
    %end
end
