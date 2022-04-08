function [Kglobal] = KmatMELAN(enrich_node,elem_crk,type_elem,xVertex,xTip,...
    split_elem,tip_elem,vertex_elem,corner_elem,tan_elem,crack_nodes,pos,xM,xCrk,Kglobal)

%declare global variables here
global node element numnode numelem elemType
global E Cm1 nu
global typeProblem typeCrack
global plotmesh plotNode
global gporder numtri
global plothelp
global orig_nn

mT = .1

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

elems = union(split_elem,vertex_elem);
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
        [flag1,flag2,crack_node] = crack_interact_element([q1,q2],e,[]);
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


  for iel=1:length(elemst)                     %loop on elems (=elements selected for enrichment)
    for kj = 1:size(xM.coor,1)-1       %loop over the elements of the fracture
      if xM.melange(kj)
        e = elemst(iel);
        q1 = xM(kk).coor(kj,:); 
        q2 = xM(kk).coor(kj+1,:);
        [flag1,flag2,crack_node] = crack_interact_element([q1,q2],e,[]);
        if flag2
            mel_elems = [mel_elems; e xM.width(kj)/2 ];
            break
        end
      end
    end
  end
end

  %loop over elements
for ii=1:size(mel_elems,1)
  % test to see if this element is within xCmel
  
  % if yes then we do all the following 
  iel = mel_elems(ii,1);
  sctr = element(iel,:) ;
  mT = mel_elems(ii,2);

  % find the distance between the two intersects (should be able to do this with det(J)
  x0 = elem_crk(iel,1) ; y0 = elem_crk(iel,2) ;
  x1 = elem_crk(iel,3) ; y1 = elem_crk(iel,4) ;
  l = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0)) ;
  nv = [(y0-y1),(x1-x0)]./l;
  mv = [(x1-x0),(y0 - y1)]./l;
  nnt = nv'*nv;
  nmt = mv'*nv; 

  JN = [ mv',nv' ];

    
  skip = 0;
  nn = length(sctr) ;
  n1 = zeros(1,nn);


  [A,BrI,QT] = f_enrich_assembly(iel,pos,type_elem,elem_crk,enrich_node);
  [A,BrI] = f_enrich_assembly(iel,pos,type_elem,elem_crk,enrich_node);
  A2 = [];
  for nI = 1:nn
    nodeI = sctr(nI);
    scB = [2*nodeI-1;2*nodeI];
    A2 = [A2;scB];
  end

  %choose Gauss quadrature rules for elements
  [W,Q] = gauss_rule(iel,enrich_node,elem_crk,...
      xTip,xVertex,tip_elem,split_elem,vertex_elem,corner_elem,xCrk) ;

  %split_corner = ismember(iel,corner_elem) & ismember(iel,split_elem)

  sctrB = [ ] ;
  for k = 1:1
      sctrB = [sctrB assembly(iel,enrich_node(:,k),pos(:,k),k,crack_nodes)] ;
  end
  scBEn = sctrB(7:end)


  %choose Gauss quadrature rules for elements
  %if strcmp(elemType,'Q4') 
    %intType = 'GAUSS' ;
  %elseif strcmp(elemType,'T3')
    %intType = 'TRIANGULAR' ;
  %end
  %[W,Q] = quadrature(2,intType,2) ;
  %loop over Gauss points
  for kk = 1:size(W,1)
      B = [] ;
      Gpt = Q(kk,:) ;
      [B, dJ] = xfemBmel(Gpt,iel,type_elem,enrich_node(:,k),elem_crk,xVertex,crack_nodes,k,JN,mT);
      % for now
      % now we want to rotate this so that is is 
      %Ppoint =  N' * node(sctr,:);
      %q = [q;Ppoint] ;

      Kglobal(A,A) = Kglobal(A,A) + B'*Cm1*B*W(kk)*dJ ;
  end
end

f = figure(1)
hold on
plotMesh(node,element(mel_elems(:,1),:),elemType,'c-','no',f)


% plot mesh with crack and enriched nodes
% plotCrack(xCrk,enrich_node,plotmesh) ;
% plot(q(:,1),q(:,2),'r*') ;
% clear q
