function [W,Q] = disVertT3(order,phi,nodes,tip,nsubDiv,intType)

global elemType
global plothelp

epsilon  = 0.000001;
corner   = [1 2 3 1];
node     = [0 0; 1 0; 0 1];


[cutEdge, node] = f_edgedetect(node, corner, phi );

vtip = f_naturalpoint(tip,nodes,20,epsilon);

% loop on element edges 
% here we have to be carefull about finding where there is an intersection with crack

% insert the tip into the Delaunay triangulation (or in this case vertex)
node = [node;vtip];

%do delaunay triangulation
tri = delaunay(node(:,1),node(:,2)) ;
tri = tricheck(node,tri) ;

if( nsubDiv == 0)
    triangles = tri ;
    triNodes = node ;
else
    [triNodes,triangles] = subTriXFEM(node,tri,nsubDiv);
end

node = triNodes ;
tri = triangles ;
% 
 %% % ---- plot of the triangulation ------------
 %% % -------------------------------------------
 %v=get(0,'ScreenSize');
 if plothelp
 %%figure('Color',[1 1 1],'Position', [0 0 0.4*v(1,3) 0.4*v(1,4)])
 figure(2);
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
    t = size(w,1) ;
%     disp(['Number of integration points     ',  num2str(t)])
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
