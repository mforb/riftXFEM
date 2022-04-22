function [ ] = f_plot_wall_forces( u,xCrk,enrDomain,typeElem,elem_force,elem_crk,split_elem,vertex_elem,tip_elem)
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




for kk = 1:size(xCrk,2) %what's the crack?
  if Hidden
    f = figure('visible','off');
  else
    f = figure();
  end
  ns = size(xCrk(kk).coor,2);
  crack_coord = {};
  f_app = {};
  for ii = 1:length(tip_elem)
    iel = tip_elem(ii);
    sctr = element(iel,:);
    segment = elem_crk(iel,:);
    for i=[1,ns-1]
      xseg = [xCrk(kk).coor(i,:),xCrk(kk).coor(i+1,:)];
      [lseg,~,~,~,~,~] = f_segment_dist(xseg);
      [flag1,flag2,~] = crack_interact_element(xseg,iel,[]);
      if flag1
        [l,~,~,~,~,~] = f_segment_dist(segment);
        dw = l/(2*sqrt(3));
        if i == 1
          crack_coord{i} = [dw, l - dw];
        else
          crack_coord{i} = [lseg - (l - dw), lseg - dw];
        end
        f_app{i} = elem_force(iel,1:2);
      end
    end
  for ii =1:length(vertex_elem)
    iel = vertex_elem(ii);
    sctr = element(iel,:);
    segment = elem_crk(iel,:);
    for i=1:ns-1
      xseg = [xCrk(kk).coor(i,:),xCrk(kk).coor(i+1,:)];
      [lseg,~,~,~,~,~] = f_segment_dist(xseg);
      [flag1,flag2,~] = crack_interact_element(xseg,iel,[]);
      if flag1
        l1 = [segment(1:2),xseg(3:4)];
        d1 = f_segment_dist(l1);
        l2 = [xseg(3:4),segment(3:4)];
        d2 = f_segment_dist(l2);
        dt = d1 + d2 ; %this is just for plotting purposes, in reality the segment between intercepts is used
        dw = dt/(2*sqrt(3));
        if d1 > dw 
          crack_coord{i} = [ crack_coord{i}, lseg - (d1 - dw) ]
          f_app{i} = [ f_app{i}, elem_force(iel,1) ]
          if d2 < dw
            crack_coord{i} = [ crack_coord{i}, lseg - d1 - (dt - dw) ]
            f_app{i} = [ f_app{i}, elem_force(iel,2) ]
          else
            crack_coord{i+1} = [ dt - dw -d1, crack_coord{i+1} ]
            f_app{i} = [ elem_force(iel,2), f_app{i+1}]
          end
        else 
          crack_coord{i+1} = [ dw-d1, dt-dw-d1 , crack_coord{i+1} ] 
          f_app{i+1} = [ elem_force(iel,2), f_pp{i+1} ] 
        end
      end
    end
  end
  for ii=1:length(split_elem)
    iel = split_elem(ii);
    sctr = element(iel,:) ;
    segment = elem_crk(iel,:);
    [l,~,~,~,~,~] = f_segment_dist(segment);
    dw = l/(2*sqrt(3));
    % if there is ocean force
    fp = elem_force(iel,1:2) ;
    for i=1:ns-1
      xseg = [xCrk(kk).coor(i,:),xCrk(kk).coor(i+1,:)];
      [flag1,flag2,~] = crack_interact_element(xseg,iel,[]);
      if flag1
        nseg = [xseg(1:2),segment(1:2)];
        [nl,~,~,~,~,~] = f_segment_dist(nseg);
        ccoords = [ nl + dw, nl + l - dw] 
        fp = elem_force(iel,:);
        seg_coord = crack_coord{i};
        seg_f = f_app{i};
        ind = find(seg_coord>ccoords(1),1);
        if isempty(ind)
          seg_coord = [seg_coord, ccoords];
          seg_f = [seg_f,fp];
        elseif ind == 1
          seg_coord = [ccoords, seg_coord];
          seg_f = [fp, seg_f];
        else
          seg_coord = [seg_coord(1:ind)];
          seg_f = [seg_f(1:ind-1),fp,seg_f(ind:end)];
        end
        crack_coord{i} = seg_coord;
        f_app{i} = seg_f;
      end
    end
  end
  pl = 0;
  cc = []; 
  ff = [];
  for i = 1:ns-1
    xseg = [xCrk(kk).coor(i,:),xCrk(kk).coor(i+1,:)];
    [lseg,~,~,~,~,~] = f_segment_dist(xseg);
    cc = [ cc, crack_coord{i} + pl ];
    ff = [ff, f_app{i}]; 
    pl = pl + lseg;
  end
  keyboard
  plot(pl,ff)
  keyboard
end




      


  end
end

