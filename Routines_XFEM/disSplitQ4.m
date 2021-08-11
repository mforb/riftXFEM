function [W,Q] = disSplitQ4(order,phi,nsubDiv,intType)

global elemType

corner = [1 2 3 4 1] ;
node = [-1 -1;1 -1;1 1;-1 1] ;

numEdge = size(node,1) ;
cutEdge = [ ] ;
%loop on element edges
for iedge = 1:numEdge
    n1 = corner(iedge) ;
    n2 = corner(iedge+1) ;
    if( phi(n1)*phi(n2) < 0 )
        r = phi(n1)/(phi(n1)-phi(n2)) ;
        pnt = (1-r)*node(n1,:)+r*node(n2,:) ;
        node = [node; pnt] ;
        %         disp(['Edge is cut by crack     ',num2str(iedge)]) ;
        cutEdge = [cutEdge iedge] ;
    end
end %end iedge loop


%check to see if adjacent edges are cut.
%If adjacent edges are, this would mean one of the sub-element will be a polygon and
%the other a triangle. Other feasibility is opposite edges are cut, which
%would mean both the sub-elements are Q4.
nEdge = length(cutEdge) ;
if(cutEdge(2) == cutEdge(1)+1 | cutEdge(2) == cutEdge(1)+3)
    %     disp('One subelement is polygon')
    if( ismember(cutEdge(1),[1 2]) & ismember(cutEdge(2),[1,2]) )   %side 1 and 2
        tempNode = [node(1,:);node(5,:);node(6,:);node(3,:);node(4,:)] ;
        [geom,iner,cpmo] = polygeom(tempNode(:,1),tempNode(:,2)) ;
        node = [node;geom(2) geom(3)] ;
        tri = [1 7 5;1 7 4;4 7 3;3 7 6;6 7 5;5 6 2] ;
        tri = tricheck(node,tri) ;
    elseif( ismember(cutEdge(1),[2 3]) & ismember(cutEdge(2),[2 3]))  %side 2 and 3
        tempNode = [node(1,:);node(2,:);node(5,:);node(6,:);node(4,:)] ;
        [geom,iner,cpmo] = polygeom(tempNode(:,1),tempNode(:,2)) ;
        node = [node;geom(2) geom(3)] ;
        tri = [1 7 2;1 7 4;4 7 6;6 7 5;5 7 2;3 6 5] ;
        tri = tricheck(node,tri) ;
    elseif( ismember(cutEdge(1),[3 4]) & ismember(cutEdge(2),[3 4]) ) %side 3 and 4
        tempNode = [node(1,:);node(2,:);node(3,:);node(5,:);node(6,:)] ;
        [geom,iner,cpmo] = polygeom(tempNode(:,1),tempNode(:,2)) ;
        node = [node;geom(2) geom(3)] ;
        tri = [1 2 7;2 7 3;7 3 5;5 6 7;6 7 1;4 5 6] ;
        tri = tricheck(node,tri) ;
    elseif( ismember(cutEdge(1),[1 4]) & ismember(cutEdge(2),[1 4]) ) %side 4 and 1
        tempNode = [node(5,:);node(2,:);node(3,:);node(4,:);node(6,:)] ;
        [geom,iner,cpmo] = polygeom(tempNode(:,1),tempNode(:,2)) ;
        node = [node;geom(2) geom(3)] ;
        tri = [1 5 6;5 7 2;6 7 5;6 7 4;4 7 3;7 3 2] ;
        tri = tricheck(node,tri) ;
    end
else
    %     disp('Both subelement are Q4')
    if( cutEdge(1) == 2 & cutEdge(2) == 4 )
        tempNode1 = [node(1,:);node(2,:);node(5,:);node(6,:)] ;
        tempNode2 = [node(6,:);node(5,:);node(3,:);node(4,:)] ;

        [geom,iner,cpmo] = polygeom(tempNode1(:,1),tempNode1(:,2)) ;
        CenPoint1 = [geom(2) geom(3)] ;

        [geom,iner,cpmo] = polygeom(tempNode2(:,1),tempNode2(:,2)) ;
        CenPoint2 = [geom(2) geom(3)] ;
        node = [node; CenPoint1; CenPoint2] ;
        tri = [1 7 2;2 5 7;7 5 6;6 7 1;6 8 5;5 8 3;3 8 4;4 8 6] ;
        tri = tricheck(node,tri) ;
    elseif( cutEdge(1) == 1 & cutEdge(2) == 3)
        tempNode1 = [node(1,:);node(5,:);node(6,:);node(4,:)] ;
        tempNode2 = [node(5,:);node(2,:);node(3,:);node(6,:)] ;

        [geom,iner,cpmo] = polygeom(tempNode1(:,1),tempNode1(:,2)) ;
        CenPoint1 = [geom(2) geom(3)] ;

        [geom,iner,cpmo] = polygeom(tempNode2(:,1),tempNode2(:,2)) ;
        CenPoint2 = [geom(2) geom(3)] ;
        node = [node; CenPoint1; CenPoint2] ;

        tri = [1 7 4;1 7 5;7 5 6;4 7 6;6 8 5;5 8 2;8 2 3;3 8 6] ;
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

% % % ---- plot of the triangulation ------------
% % % -------------------------------------------
% v=get(0,'ScreenSize');
% %figure('Color',[1 1 1],'Position', [0 0 0.4*v(1,3) 0.4*v(1,4)])
% nd=[];
%    for igp = 1 : size(node,1)
%         gpnt = node(igp,:);
%         [N,dNdxi]=lagrange_basis(elemType,gpnt);
%         Gpnt = N' * nodes; % global GP
%         nd = [nd;Gpnt];
%    end
% triplot(tri, nd(:,1),nd(:,2),'g')
% 
% % ------------------------------------------

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
