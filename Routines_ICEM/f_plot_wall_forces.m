function [ ] = f_plot_wall_forces( u,xCrk,enrDomain,typeElem,elem_force,elem_gap,elem_crk,split_elem,vertex_elem,tip_elem,stepnum)
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
  crack_gapc = {};
  f_app = {};
  f_tra = {};
  gapn  = {};
  gapt  = {};
  gapw  = {};
  if ns == 2
    f_app{1} = [0,0] ;
    f_tra{1} = [0,0] ; 
    gapn{1}  = [0,0] ;
    gapt{1}  = [0,0] ;
    gapw{1}  = [0,0] ;
    crack_coord{1} = [0 lseg]; 
    crack_gapc{1} = [0 lseg]; 
  else
    f_app{1} = 0; f_app{ns-1} = 0;
    f_tra{1} = 0; f_tra{ns-1} = 0;
    gapn{1}  = 0; gapn{ns-1} = 0;
    gapt{1}  = 0; gapt{ns-1} = 0;
    gapw{1}  = 0 ; gapw{ns-1} = 0;
    crack_coord{1} = 0; crack_coord{ns-1} = lseg;
    crack_gapc{1} = 0; crack_gapc{ns-1} = lseg;
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
            crack_gapc{i} = [crack_gapc{i},l];
            gapn{i}  = [ gapn{i},elem_gap(iel,4)];
            gapt{i}  = [ gapt{i},elem_gap(iel,3)];
          else
            ccoord = lseg-ccoord; 
            [ccoord,so] = sort(ccoord);
            crack_coord{i} = [ccoord, crack_coord{i}];
            crack_gapc{i} = [lseg-l,crack_gapc{i}];
            f_app{i} = [  fp(so), f_app{i}];
            f_tra{i} = [  ft(so), f_tra{i}];
            gapn{i}  = [ elem_gap(iel,2), gapn{i}];
            gapt{i}  = [ elem_gap(iel,1), gapt{i}];
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
    [~,width] = f_find_melange(iel,xCrk(kk));
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
        gc     = [nl,nl+l]
        fp = elem_force(1,iel,1:2:wall_int*2); fp = fp(so); fp = reshape(fp,1,wall_int);
        ft = elem_force(1,iel,2:2:wall_int*2); ft = ft(so); ft = reshape(ft,1,wall_int);
        %fp = (fp)/l;
        %ft = (ft)/l;
        seg_coord = crack_coord{i};
        seg_gc    = crack_gapc{i};
        seg_f = f_app{i};
        seg_t = f_tra{i};
        seg_gn = gapn{i};
        seg_gt = gapt{i};
        seg_w = gapw{i};
        ind = find(seg_coord>ccoord(end),1);
        if isempty(ind)
          seg_coord = [seg_coord, ccoord];
          seg_gc   = [gc,seg_gc];
          seg_f = [seg_f,fp];
          seg_t = [seg_t,ft];
          seg_gn = [seg_gn,elem_gap(iel,[2,4])];
          seg_gt = [seg_gt,elem_gap(iel,[1,3])];
          seg_w = [seg_w,width,width];
        elseif ind == 1
          seg_coord = [ccoord, seg_coord];
          seg_gc   = [seg_gc,gc];
          seg_f = [seg_f,fp];
          seg_f = [fp, seg_f];
          seg_t = [ft, seg_t];
          seg_gn = [elem_gap(iel,[2,4]),seg_gn];
          seg_gt = [elem_gap(iel,[1,3]),seg_gt];
          seg_w = [width,width,seg_w];
        elseif ind == 1
        else
          seg_coord = [seg_coord(1:ind-1),ccoord,seg_coord(ind:end)];
          seg_gc   = [seg_gc(1:ind-1),gc,seg_gc(ind:end)];
          seg_f = [seg_f(1:ind-1),fp,seg_f(ind:end)];
          seg_t = [seg_t(1:ind-1),ft,seg_t(ind:end)];
          seg_gn = [seg_gn(1:ind-1),elem_gap(iel,[2,4]),seg_gn(ind:end)];
          seg_gt = [seg_gt(1:ind-1),elem_gap(iel,[1,3]),seg_gt(ind:end)];
          seg_w = [seg_w(1:ind-1),width,width,seg_w(ind:end)];
        end
        crack_coord{i} = seg_coord;
        f_app{i} = seg_f;
        f_tra{i} = seg_t;
        crack_gapc{i} = seg_gc;
        gapn{i} = seg_gn;
        gapt{i} = seg_gt;
        gapw{i} = seg_w;
      end
    end
  end
  for ii =1:length(vertex_elem)
    iel = vertex_elem(ii);
    sctr = element(iel,:);
    segment = elem_crk(iel,:);
    [~,width] = f_find_melange(iel,xCrk(kk));
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
        gc1 = [lseg-l1,lseg];
        gc2 = [0,l2];
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
        crack_gapc{i} = [ crack_gapc{i}, gc1 ];
        gapn{i} = [ gapn{i}, elem_gap(iel,[2,4]) ];
        gapt{i} = [ gapt{i}, elem_gap(iel,[1,3]) ];
        gapw{i} = [ gapw{i}, width,width ];
        crack_coord{i+1} = [ ccoord2, crack_coord{i+1} ];
        f_app{i+1} = [ fp, f_app{i+1}];
        f_tra{i+1} = [ ft, f_tra{i+1}];
        crack_gapc{i+1} = [ gc2, crack_gapc{i+1} ];
        gapn{i+1} = [ elem_gap(iel,[4,6]), gapn{i+1} ];
        gapt{i+1} = [ elem_gap(iel,[3,5]),  gapt{i+1} ];
        gapw{i+1} = [ elem_gap(iel,[3,5]), gapw{i+1}  ];
        break
      end
    end
  end
  pl = 0;
  cc = []; 
  gc = [];
  ff = [];
  tt = [];
  gn = [];
  gt = [];
  gw = [];
  inters = [];
  for i = 1:ns-1
    xseg = [xCrk(kk).coor(i,:),xCrk(kk).coor(i+1,:)];
    [lseg,~,~,~,~,~] = f_segment_dist(xseg);
    cc = [ cc, crack_coord{i} + pl ];
    gc = [ gc, crack_gapc{i} + pl ];
    ff = [ff, f_app{i}]; 
    tt = [tt, f_tra{i}]; 
    gn = [gn, gapn{i}]; 
    gt = [gt, gapt{i}]; 
    gw = [gw, gapw{i}]; 
    pl = pl + lseg;
    inters = [inters, pl];
  end
  inters(end) = [];

  t = tiledlayout(2,1,'TileSpacing','Compact');
  nexttile

  plot(cc,ff,'linewidth',3,'DisplayName','normal sym force')
  yl = ylim();
  hold on
  for i = 1: length(inters)
    plot([inters(i),inters(i)],yl,'c-','linewidth',0.8,'Color',[0.5,0.3,0.7,0.2]);
  end
  ylim(yl);
  xlabel('distance along crack')
  ylabel('force')
  tstr = ['normal forces acting on rift ',num2str(kk)];
  title(tstr);

  nexttile
  plot(cc,tt,'color',[1,0.1,0],'linewidth',2,'DisplayName','tangential force')
  hold on
  yl = ylim();
  for i = 1: length(inters)
    plot([inters(i),inters(i)],yl,'c-','linewidth',0.8,'Color',[0.5,0.3,0.7,0.2]);
  end
  ylim(yl);
  xlabel('distance along crack')
  ylabel('force')
  tstr = ['tangential forces acting on rift ',num2str(kk)];
  title(tstr);
  nstr = ['rift',num2str(kk),'forces_step',num2str(stepnum)];
  print([results_path,'/',nstr],'-dpng','-r300')
  clf


  t = tiledlayout(3,1,'TileSpacing','Compact');
  nexttile

  plot(gc,gn,'linewidth',3,'DisplayName','normal gap')
  yl = ylim();
  hold on
  for i = 1: length(inters)
    plot([inters(i),inters(i)],yl,'c-','linewidth',0.8,'Color',[0.5,0.3,0.7,0.2]);
  end
  ylim(yl);
  xlabel('distance along crack (m)')
  ylabel('gap (m)')
  tstr = ['normal gap along rift ',num2str(kk)];
  title(tstr);

  nexttile
  plot(gc,gt,'color',[1,0.1,0],'linewidth',2,'DisplayName','tangential disp')
  hold on
  yl = ylim();
  for i = 1: length(inters)
    plot([inters(i),inters(i)],yl,'c-','linewidth',0.8,'Color',[0.5,0.3,0.7,0.2]);
  end
  ylim(yl);
  xlabel('distance along crack (m)')
  ylabel('slip (m)')
  tstr = ['tangential displacement along rift ',num2str(kk)];
  title(tstr);

  nexttile
  plot(gw+gn,gt,'color',[0,0.7,0.1],'linewidth',2,'DisplayName','rift wall distance')
  hold on
  yl = ylim();
  for i = 1: length(inters)
    plot([inters(i),inters(i)],yl,'c-','linewidth',0.8,'Color',[0.5,0.3,0.7,0.2]);
  end
  ylim(yl);
  xlabel('distance along crack (m)')
  ylabel('distance (m)')
  tstr = ['"melange" distance between rift walls ',num2str(kk)];
  title(tstr);



  nstr = ['rift',num2str(kk),'gaps_step',num2str(stepnum)];
  print([results_path,'/',nstr],'-dpng','-r300')
end

