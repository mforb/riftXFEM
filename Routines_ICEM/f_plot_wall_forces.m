function [ ] = f_plot_wall_forces( u,xCrk,enrDomain,typeElem,elem_force,elem_crk,split_elem,vertex_elem,tip_elem,stepnum)
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Fri Mar 18 18:39:27 NZDT 2022

%declare global variables here
global node element numnode numelem elemType
global plotmesh plotNode
global gporder numtri
global plothelp Hidden results_path
global rift_wall_pressure
global wall_int



%loop over elements
%elems = union(split_elem,vertex_elem);


[W,Q] = quadrature(wall_int,'GAUSS',1) ;

for kk = 1:size(xCrk,2) %what's the crack?
  if Hidden
    f = figure('visible','off');
  else
    f = figure();
  end
  ns = size(xCrk(kk).coor,1);
  xseg = [xCrk(kk).coor(end-1,:),xCrk(kk).coor(end,:)];
  [lseg,~,~,~,~,~] = f_segment_dist(xseg);
  crack_coord = {};
  f_app = {};
  f_tra = {};
  if ns == 2
    f_app{1} = [0,0] ;
    f_tra{1} = [0,0] ; 
    crack_coord{1} = [0 lseg]; 
  else
    f_app{1} = 0; f_app{ns-1} = 0;
    f_tra{1} = 0; f_tra{ns-1} = 0;
    crack_coord{1} = 0; crack_coord{ns-1} = lseg;
  end
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
        ccoord = ((l/2).*(Q+1))';
        [ccoord,so] = sort(ccoord);
        fp = elem_force(1,iel,1:2:wall_int*2); fp=fp(so); fp = reshape(fp,1,wall_int);
        ft = elem_force(1,iel,2:2:wall_int*2); ft = ft(so); ft = reshape(ft,1,wall_int);
        %fp = (fp)/l;
        %ft = (ft)/l;
        if ns == 2
          if ~points_same_2d(segment(1:2),xseg(1:2))
            ccoord = lseg-ccoord; 
          end
          [ccoord,so] = sort(ccoord);
          ind = find(crack_coord{1}>ccoord(end),1);
          crack_coord{1}= [crack_coord{1}(1:ind-1),ccoord,crack_coord{1}(ind:end)];
          f_app{1} = [f_app{1}(1:ind-1),fp(so),f_app{1}(ind:end)];
          f_tra{1} = [f_tra{1}(1:ind-1),ft(so),f_tra{1}(ind:end)];
          break
        else
          if i == 1
            crack_coord{i} = [crack_coord{i}, ccoord];
            f_app{i} = [ f_app{i}, fp ];
            f_tra{i} = [ f_tra{i}, ft ];
          else
            ccoord = lseg-ccoord; 
            [ccoord,so] = sort(ccoord);
            crack_coord{i} = [ccoord, crack_coord{i}];
            f_app{i} = [  fp(so), f_app{i}];
            f_tra{i} = [  ft(so), f_tra{i}];
          end
        end
      end
    end
  end
  [~,so] = sort(Q);
  for ii=1:length(split_elem)
    iel = split_elem(ii);
    sctr = element(iel,:) ;
    segment = elem_crk(iel,:);
    [l,~,~,~,~,~] = f_segment_dist(segment);
    % if there is ocean force
    for i=1:ns-1
      xseg = [xCrk(kk).coor(i,:),xCrk(kk).coor(i+1,:)];
      [flag1,flag2,~] = crack_interact_element(xseg,iel,[]);
      [~,flag_int] = f_edge_int_points(iel,xseg);
      if flag1 & flag_int
        nseg = [xseg(1:2),segment(1:2)];
        [nl,~,~,~,~,~] = f_segment_dist(nseg);
        ccoord = (l/2).*(Q+1)'+nl;
        ccoord = ccoord(so);
        fp = elem_force(1,iel,1:2:wall_int*2); fp = fp(so); fp = reshape(fp,1,wall_int);
        ft = elem_force(1,iel,2:2:wall_int*2); ft = ft(so); ft = reshape(ft,1,wall_int);
        %fp = (fp)/l;
        %ft = (ft)/l;
        seg_coord = crack_coord{i};
        seg_f = f_app{i};
        seg_t = f_tra{i};
        ind = find(seg_coord>ccoord(end),1);
        if isempty(ind)
          seg_coord = [seg_coord, ccoord];
          seg_f = [seg_f,fp];
          seg_t = [seg_t,ft];
        elseif ind == 1
          seg_coord = [ccoord, seg_coord];
          seg_f = [fp, seg_f];
          seg_t = [ft, seg_t];
        else
          seg_coord = [seg_coord(1:ind-1),ccoord,seg_coord(ind:end)];
          seg_f = [seg_f(1:ind-1),fp,seg_f(ind:end)];
          seg_t = [seg_t(1:ind-1),ft,seg_t(ind:end)];
        end
        crack_coord{i} = seg_coord;
        f_app{i} = seg_f;
        f_tra{i} = seg_t;
      end
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
      [~,flag_int] = f_edge_int_points(iel,xseg);
      if flag1 & flag_int
        nseg = [xseg(1:2),segment(1:2)];
        [nl,~,~,~,~,~] = f_segment_dist(nseg);
        l1 = [segment(1:2),xseg(3:4)];
        d1 = f_segment_dist(l1);
        l2 = [xseg(3:4),segment(3:4)];
        d2 = f_segment_dist(l2);
        ccoord1 = (d1/2).*(Q+1)'+nl;
        ccoord1 =ccoord1(so);
        ccoord2 = (d2/2).*(Q+1)';
        ccoord2 =ccoord2(so);
        fp1 = elem_force(1,iel,1:2:wall_int*2); fp1 = fp1(so); fp1 = reshape(fp1,1,wall_int);
        ft1 = elem_force(1,iel,2:2:wall_int*2); ft1 = ft1(so); ft1 = reshape(ft1,1,wall_int);
        %fp1 = (fp1)/d1;
        %ft1 = (ft1)/d1;
        fp2 = elem_force(2,iel,1:2:wall_int*2); fp2 = fp2(so); fp2 = reshape(fp2,1,wall_int);
        ft2 = elem_force(2,iel,2:2:wall_int*2); ft2 = ft2(so); ft2 = reshape(ft2,1,wall_int);
        %fp2 = (fp2)/d1;
        %ft2 = (ft2)/d2;
        crack_coord{i} = [ crack_coord{i}, ccoord1 ];
        f_app{i} = [ f_app{i}, fp1 ];
        f_tra{i} = [ f_tra{i}, ft1 ];
        crack_coord{i+1} = [ ccoord2, crack_coord{i+1} ];
        f_app{i+1} = [ fp, f_app{i+1}];
        f_tra{i+1} = [ ft, f_tra{i+1}];
        break
      end
    end
  end
  pl = 0;
  cc = []; 
  ff = [];
  tt = [];
  for i = 1:ns-1
    xseg = [xCrk(kk).coor(i,:),xCrk(kk).coor(i+1,:)];
    [lseg,~,~,~,~,~] = f_segment_dist(xseg);
    cc = [ cc, crack_coord{i} + pl ];
    ff = [ff, f_app{i}]; 
    tt = [tt, f_tra{i}]; 
    pl = pl + lseg;
  end

  t = tiledlayout(2,1,'TileSpacing','Compact');
  nexttile

  plot(cc,ff,'linewidth',3,'DisplayName','normal sym force')
  xlabel('distance along crack')
  ylabel('force')
  tstr = ['normal forces acting on rift ',num2str(kk)];
  title(tstr);

  nexttile
  plot(cc,tt,'color',[1,0.1,0],'linewidth',2,'DisplayName','tangential force')
  xlabel('distance along crack')
  ylabel('force')
  tstr = ['tangential forces acting on rift ',num2str(kk)];
  title(tstr);
  nstr = ['rift',num2str(kk),'forces_step',num2str(stepnum)];
  print([results_path,'/',nstr],'-dpng','-r300')
  clf
end

