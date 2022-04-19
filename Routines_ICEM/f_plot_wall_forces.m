function [ ] = f_plot_wall_forces( u,xCr,enrDomain,typeElem,elem_crk,split_elem)
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Fri Mar 18 18:39:27 NZDT 2022

%declare global variables here
global node element numnode numelem elemType
global plotmesh plotNode
global gporder numtri
global plothelp
global rift_wall_pressure



%loop over elements
elems = union(split_elem,vertex_elem);

for kk = 1:size(xCrk,2) %what's the crack?
  for ii=1:size(elems,1)
    iel = elems(ii);
    sctr = element(iel,:) ;
    segment = elem_crk(iel,:);
    x0 = elem_crk(iel,1) ; y0 = elem_crk(iel,2) ;
    x1 = elem_crk(iel,3) ; y1 = elem_crk(iel,4) ;
    l = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0)) ;
    nv = [(y0-y1),(x1-x0)]./l;
    % if there is ocean force
    f_app = elem_force(iel,1:2) ;

    % need to sort the segments into one coordinate (along crack)
    seg_value = 0
    coord = 0
    seg_cumul = 0
    for 
      [flag1,flag2,~] = crack_interact_element(q,e,[])


  end
end

