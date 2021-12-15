function [W,Q] = disSplitT3(order,phi,nodes,nsubDiv,intType,hascorner)

global elemType
global plothelp
global orig_nn

corner = [1 2 3 1] ;
node = [0 0;1 0;0 1] ;

[cutEdge, node] = f_edgedetect(node, corner,  phi) ;

%check to see if adjacent edges are cut.
%If adjacent edges are, this would mean one of the sub-element will be a polygon and
%the other a triangle. Other option is that a corner is split this means that both sides are triangles  
nEdge = length(cutEdge) ;


if ~hascorner
  if( ismember(cutEdge(1),[1 2]) & ismember(cutEdge(2),[1,2]) )   %side 1 and 2
      tempNode = [node(1,:);node(4,:);node(5,:);node(3,:)] ;
      [geom,iner,cpmo] = polygeom(tempNode(:,1),tempNode(:,2)) ;
      node = [node;geom(2) geom(3)] ;
      tri = [1 6 4;1 6 3;4 6 5;3 6 5;4 5 2] ;
      tri = tricheck(node,tri) ;
  elseif( ismember(cutEdge(1),[2 3]) & ismember(cutEdge(2),[2 3]))  %side 2 and 3
      tempNode = [node(1,:);node(2,:);node(4,:);node(5,:)] ;
      [geom,iner,cpmo] = polygeom(tempNode(:,1),tempNode(:,2)) ;
      node = [node;geom(2) geom(3)] ;
      tri = [1 6 2;1 6 5;5 6 4;2 6 4;3 5 4] ;
      tri = tricheck(node,tri) ;
  elseif( ismember(cutEdge(1),[3 1]) & ismember(cutEdge(2),[3 1]) ) %side 3 and 1
      tempNode = [node(4,:);node(2,:);node(3,:);node(5,:)] ;
      [geom,iner,cpmo] = polygeom(tempNode(:,1),tempNode(:,2)) ;
      node = [node;geom(2) geom(3)] ;
      tri = [2 4 6;2 6 3;4 6 5;5 6 3;1 4 5] ;
      tri = tricheck(node,tri) ;
  end
else
  if cutEdge == 1
    tri = [ 1 4 3; 2 3 4 ];
    tri = tricheck(node,tri) ;
  elseif cutEdge == 2
    tri = [ 1 4 3; 2 4 1 ];
    tri = tricheck(node,tri) ;
  elseif cutEdge == 3
    tri = [ 1 2 4; 2 3 4 ];  
    tri = tricheck(node,tri) ;
  end
end



if( nsubDiv == 0)
    triangles = tri ;
    triNodes = node ;
else
    [triNodes,triangles] = subTriXFEM(node,tri,nsubDiv);
end

node = triNodes ;
tri = triangles ;

 %% % ---- plot of the triangulation ------------
 %% % -------------------------------------------
 %v=get(0,'ScreenSize');
 %%figure('Color',[1 1 1],'Position', [0 0 0.4*v(1,3) 0.4*v(1,4)])
if plothelp
figure(2);
hold on
nd=[];
    for igp = 1 : size(node,1)
         gpnt = node(igp,:);
         [N,dNdxi]=lagrange_basis(elemType,gpnt);
         Gpnt = N' * nodes; % global GP
         nd = [nd;Gpnt];
    end
triplot(tri, nd(:,1),nd(:,2),'g--')
end
 
 %% ------------------------------------------

% loop over subtriangles to get quadrature points and weights
pt = 1;
for e=1:size(tri,1)
    [w,q] = quadrature(order,intType,2);
    % transform quadrature points into the parent element
    coord = node(tri(e,:),:);
    a = det([coord,[1;1;1]])/2;
    if ( a<0 )  % need to swap connectivity
        coord=[coord(2,:);coord(1,:);coord(3,:)];
        a = det([coord,[1;1;1]])/2;
    end

    if ( a~=0 )
        for n=1:length(w)
            N = lagrange_basis('T3',q(n,:));
            Q(pt,:) = N'*coord;
            W(pt,1) = 2*w(n)*a;
            pt = pt+1;
        end
    end

end
