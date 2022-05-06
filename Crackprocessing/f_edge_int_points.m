function [int_points,flag_int] = f_edgede_int_points(iel,xseg);
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Mon May 2 2022

global element node

sctr=element(iel,:);
sctrl = [sctr sctr(1,1)];      
int_points = [];
flag_int = 0;
for iedge=1:size(sctr,2)             %loop over the edges of elements
    nnode1=sctrl(iedge);
    nnode2=sctrl(iedge+1);
    p1 = [node(nnode1,:)];
    p2 = [node(nnode2,:)];
    intersect = segments_int_2d(p1,p2,xseg(1:2),xseg(3:4)) ;
    %edge_track=[edge_track,iedge]
    if intersect(1)                
      flag_int = 1;
      int_points = [ int_points ; intersect(2:3) ];
    end
end
