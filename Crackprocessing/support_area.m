function [Aw, Awp] = support_area(pt,e,type_elem,elem_crk,xVertex,kk)

% function made for all nodes
%
% pt = is the number of the point 
% e = number of THE considered element (see node detect) 
% type_elem = vector which give the type of all elements (1=tip, 2=split,
%                                               3=vertex,0=normal)
% xCr1 = coordinate of all the crack intersection with the mesh
% kk = index of the crack (for multiple crack)  


global node element

[sctrn,xx] = find(element == pt);                   % find support elements
coor_pt = node(pt,:);                               %
xCre = [elem_crk(e,1) elem_crk(e,2) ;...
        elem_crk(e,3) elem_crk(e,4)];               % intersection with the element
distPT = signed_distance(xCre,coor_pt,1) ;          % below or above?  (sign)
HPT = sign(distPT)  ;                               % below or above?  (-1 or 1)




Awp = 0; % area +
Aw  = 0; % total area

for i = 1 : size(sctrn,1)                                    
    sctr = element(i,:);                                     
    Aw = Aw + polyarea(node(sctr,1),node(sctr,2));            
end               

if (HPT == 0) | (abs(distPT) < 1e-4)                % la fessura coincide col bnordo di un elemento!
    HPT = 1;
%     disp('the node coincide with the crack')
%     return
end

for i = 1: size(sctrn,1)                            % loop on the elements of the support
     sctr = element(sctrn(i),:);
     nn   = length(sctr);
     pts = [];
     
     if (type_elem(sctrn(i),kk) == 0)
         for in = 1:nn                              % loop on the nodes of the element
             nd = node(sctr(in),:);
             dist = signed_distance(xCre,nd,1);
             Hi = sign(dist);
             if Hi == HPT                           % is the in th same side of the node pt
                 pts = [pts; nd];
             end
         end
     end
 
     if (type_elem(sctrn(i),kk) == 2) | (type_elem(sctrn(i),kk) == 3)               % is a cut element 
        xCr_local = [elem_crk(sctrn(i),1) elem_crk(sctrn(i),2);...
                      elem_crk(sctrn(i),3) elem_crk(sctrn(i),4)];   
        for in = 1:nn                                   % dist again to be sure of the side of each point
            nd = node(sctr(in),:);
            dist = signed_distance(xCr_local,nd,1);
            Hi = sign(dist);
            if Hi == HPT                                % is the in th same side of the node pt
                pts = [pts; nd];
            end
        end
        pts = [pts; xCr_local] ;                              % now all points of the polygonal Awp domain
    end
 
    if (size(pts,1) > 2 )                               % test : is this a polygonal element
        [xx,surf] = convhull(pts(:,1),pts(:,2));        % compute the surface of the convexe area
        Awp = Awp + surf;
    else    
        Awp = Awp;
    end
    
    if (type_elem(sctrn(i),kk) == 3)                    % concerning the triangle formed by the crack ends and the vertex
        dist = signed_distance(xCr_local,xVertex(sctrn(i),:),1);
        Hi = sign(dist);
        points = [xCr_local ; xVertex(sctrn(i),:)];
        [xx,surf] = convhull(points(:,1),points(:,2));
        Awp = Awp - Hi*HPT*surf ;                       % add or substract the area of this triangle
    end       
end

