function [Kglobal,nodeTanfix] = KmatMELAN(enrich_node,elem_crk,type_elem,xVertex,xTip,...
    split_elem,tip_elem,vertex_elem,corner_elem,tan_elem,crack_node,pos,xM,xCrk,Kglobal,nodeTanfix)

%declare global variables here
global node element numnode numelem elemType
global E Cm1 nu
global typeProblem typeCrack
global plotmesh plotNode
global gporder numtri
global plothelp
global orig_nn


if strcmp(elemType,'Q4') 
  intType = 'GAUSS' ;
elseif strcmp(elemType,'T3')
  intType = 'TRIANGULAR' ;
end

elems = union(split_elem,vertex_elem);
elemst = tan_elem;

mel_elems = []
for kk = 1:size(xCrk,2)
  for iel=1:length(elems)                     %loop on elems (=elements selected for enrichment)
    [flag1,width,nodeTanfix] = f_find_melange(elems(iel),xCrk,nodeTanfix);
    if flag1
      mel_elems = [mel_elems; elems(iel), width ];
    end
  end
end

keyboard


  %loop over elements
for ii=1:size(mel_elems,1)
  % test to see if this element is within xCmel
  
  % if yes then we do all the following 
  iel = mel_elems(ii,1);
  sctr = element(iel,:) ;
  nn = length(sctr) ;
  mT = mel_elems(ii,2);

  % find the distance between the two intersects (should be able to do this with det(J)
  [l,nv,mv,nnt,nmt,mmt] = f_segment_dist(elem_crk(iel,:));

  JN = [ mv',nv' ]; % this is the rotation matrix
  [A,~,~,~,~] = f_enrich_assembly(iel,pos,type_elem,elem_crk,enr_node);
  [W,Q] = quadrature(2,intType,2) ;



  for kk = 1:size(W,1)
      B = [] ;
      Gpt = Q(kk,:) ;
      [B, dJ] = xfemBmel(Gpt,iel,type_elem,enrich_node(:,k),elem_crk,xVertex,crack_node,k,JN,mT);
      % for now
      % now we want to rotate this so that is is 
      %Ppoint =  N' * node(sctr,:);
      %q = [q;Ppoint] ;

      Kglobal(A,A) = Kglobal(A,A) + B'*Cm1*B*W(kk)*dJ ;
  end
end
