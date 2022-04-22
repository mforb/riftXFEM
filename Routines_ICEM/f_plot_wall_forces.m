function [ ] = f_plot_wall_forces( u,xCr,enrDomain,typeElem,elem_crk,split_elem,vertex_elem,tip_elem)
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Fri Mar 18 18:39:27 NZDT 2022

%declare global variables here
global node element numnode numelem elemType
global plotmesh plotNode
global gporder numtri
global plothelp Hidden
global rift_wall_pressure



%loop over elements
%elems = union(split_elem,vertex_elem);
if 
  f = figure('visible','off');
else
  f = figure();
end




for kk = 1:size(xCrk,2) %what's the crack?
  if 
    f = figure('visible','off');
  else
    f = figure();
  end
  ns = size(xCrk(kk).coor,2);
  crack_coord = {};
  f_app = {};
  for ii = 1:length(tip_elems)
    iel = tip_elem(ii);
    sctr = element(iel,:);
    segment = elem_crk(iel,:);
    for i=[1,ns-1]
      xseg = [xCrk(kk).coor(i,:),xCrk(kk).coor(i+1,:)];
      [lseg,~,~,~,~,~] = f_segment_dist(xseg);
      [flag1,flag2,~] = crack_interact_element(xseg,e,[])
      if flag1
        [l,~,~,~,~,~] = f_segment_dist(seg);
        if i == 1
          crack_coord{i} = [l/(2*sqrt(3), l - l/(2*sqrt(3)];
        else
          crack_coord{i} = [lseg - (l -l/(2*sqrt(3)), lseg - l/(2*sqrt(3)];
        end
        f_app{i} = elem_forces(iel,1:2);
      end
    end
  for ii =1:length(vertex_elem)
    iel = vertex_elem(ii);
    sctr = element(iel,:);
    segment = elem_crk(iel,:);
    for i=1:ns-1
      xseg = [xCrk(kk).coor(i,:),xCrk(kk).coor(i+1,:)];
      [lseg,~,~,~,~,~] = f_segment_dist(xseg);
      [flag1,flag2,~] = crack_interact_element(xseg,e,[])
      if flag1
        l1 = [segment(1:2),xseg(3:4)];
        d1 = f_segment_dist(l1);
        l2 = [xseg(3:4),segment(3:4)];
        d2 = f_segment_dist(l2);
        dt = d1 + d2 ; %this is just for plotting purposes, in reality the segment between intercepts is used
        dw = dt/(2*sqrt(3));
        if d1 > dw 
          crack_coord{i} = [ crack_coord{i}, lseg - (d1 - dw) ]
          f_app{i} = [ f_app{i}, elem_forces(iel,1) ]
          if d2 < dw
            crack_coord{i} = [ crack_coord{i}, lseg - d1 - (dt - dw) ]
            f_app{i} = [ f_app{i}, elem_forces(iel,2) ]
          else
            crack_coord{i+1} = [ dt - dw -d1, crack_coord{i+1} ]
            f_app{i} = [ elem_forces(iel,2), f_app{i+1}]
        else 
          crack_coord{i+1} = [ dw-d1, dt-dw-d1 , crack_coord{i+1} ] 
          f_app{i+1} = [ elem_forces(iel,2), f_pp{i+1} ] 
        end
      end
    end
  end
  for ii=1:length(split_elems)
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
    for i=1:(size(xCrk(kk).coor,2)-1)
      xseg = [xCrk(kk).coor(i,:),xCrk(kk).coor(i+1,:)];
    for 
      [flag1,flag2,~] = crack_interact_element(xseg,e,[])
      


  end
end

