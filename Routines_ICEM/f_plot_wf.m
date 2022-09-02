function [ ] = f_plot_wf( u,xCrk,enrDomain,typeElem,elem_force,elem_gap,elem_crk,split_elem,vertex_elem,tip_elem,stepnum)
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Fri Mar 18 18:39:27 NZDT 2022

%declare global variables here
global node element numnode numelem elemType
global plotmesh plotNode
global gporder numtri
global plothelp Hidden results_path
global rift_wall_pressure
global wall_int
global melange




%loop over elements
%elems = union(split_elem,vertex_elem);

if isempty(wall_int)
  wall_int = 2;
end

if isempty(melange)  | melange == 0 
  mel_b = 0;
else
  mel_b = 1
end



[W,Q] = quadrature(wall_int,'GAUSS',1) ;

for kk = 1:size(xCrk,2) %what's the crack?
  if Hidden
    f = figure('visible','off');
  else
    f = figure();
  end
  ns = size(xCrk(kk).coor,1);
  if ns == 2
    ids = [1];
  else
    ids = [1,ns-1];
  end
  xseg = [xCrk(kk).coor(end-1,:),xCrk(kk).coor(end,:)];
  [lseg,~,~,~,~,~] = f_segment_dist(xseg);
  crack_coord = {[]};
  crack_gapc = {[]};
  f_app = {[]};
  f_tra = {[]};
  gapn  = {[]};
  gapt  = {[]};
  gapw  = {[]};
  for ii = 1:length(tip_elem)
    iel = tip_elem(ii);
    sctr = element(iel,:);
    segment = elem_crk(iel,:);

    if mel_b
      [~,width] = f_find_melange(iel,xCrk(kk));
    else
      width = 0;
    end
    
    for i=ids

      if i~=1
        crack_coord{i} = [];
        crack_gapc{i} = [];
        f_app{i} = [];
        f_tra{i} = [];
        gapn{i} = [];
        gapt{i} = [];
        gapw{i}  = [];
      end

      xseg = [xCrk(kk).coor(i,:),xCrk(kk).coor(i+1,:)];
      [lseg,~,~,~,~,~] = f_segment_dist(xseg);
      [flag1,flag2,~] = crack_interact_element(xseg,iel,[]);
      if flag1
        [l,~,~,~,~,~] = f_segment_dist(segment);
        ccoord = ((l/2).*(Q+1))';
        [ccoord,so] = sort(ccoord);
        ccoord = [0,ccoord]; % for the ends
        gcoord = [0,l];
        ind = 4;
        fp = elem_force(1,iel,1:2:wall_int*2); fp=fp(so); fp = reshape(fp,1,wall_int);
        ft = elem_force(1,iel,2:2:wall_int*2); ft = ft(so); ft = reshape(ft,1,wall_int);
        if ~points_same_2d(segment(1:2),xseg(1:2))
          ccoord = lseg - l + ccoord; 
          ccoord(1) = lseg; % the first coordinate is still the end
          gcoord = lseg - gcoord;
          ind = 2;
        end
        crack_coord{i} = [crack_coord{i}, ccoord];
        f_app{i} = [ f_app{i},0,fp ];
        f_tra{i} = [ f_tra{i},0,ft ];
        crack_gapc{i} = [crack_gapc{i},gcoord];
        gapn{i}  = [ gapn{i},0,elem_gap(iel,ind)];
        gapt{i}  = [ gapt{i},0,elem_gap(iel,ind-1)];
        gapw{i}  = [ gapw{i},0,2*width];
      end
    end
  end
  % Now we do all the split elements
  for ii=1:length(split_elem)
    iel = split_elem(ii);
    sctr = element(iel,:) ;
    segment = elem_crk(iel,:);
    [l,~,~,~,~,~] = f_segment_dist(segment);
    if mel_b
      [~,width] = f_find_melange(iel,xCrk(kk));
    else
      width = 0;
    end

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
        gc     = [nl,nl+l];
        fp = elem_force(1,iel,1:2:wall_int*2); fp = fp(so); fp = reshape(fp,1,wall_int);
        ft = elem_force(1,iel,2:2:wall_int*2); ft = ft(so); ft = reshape(ft,1,wall_int);
        crack_coord{i} = [crack_coord{i}, ccoord];
        crack_gapc{i}   = [crack_gapc{i},gc];
        f_app{i} =  [f_app{i},fp];
        f_tra{i} = [f_tra{i},ft];
        gapn{i} = [gapn{i},elem_gap(iel,[2,4])];
        gapt{i} = [gapt{i},elem_gap(iel,[1,3])];
        gapw{i} = [gapw{i},width,width];
      end
    end
  end
  % now the vertex elements
  for ii =1:length(vertex_elem)
    iel = vertex_elem(ii);
    sctr = element(iel,:);
    segment = elem_crk(iel,:);
    if mel_b
      [~,width] = f_find_melange(iel,xCrk(kk));
    else
      width = 0;
    end
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
        gc1 = [lseg-d1,lseg];
        gc2 = [0,d2];
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
        f_app{i+1} = [ fp2, f_app{i+1}];
        f_tra{i+1} = [ ft2, f_tra{i+1}];
        crack_gapc{i+1} = [ gc2, crack_gapc{i+1} ];
        gapn{i+1} = [ elem_gap(iel,[4,6]), gapn{i+1} ];
        gapt{i+1} = [ elem_gap(iel,[3,5]),  gapt{i+1} ];
        gapw{i+1} = [ width,width, gapw{i+1}  ];
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
    [tcc,so1] = sort(crack_coord{i});
    [tgc,so2] = sort(crack_gapc{i});
    cc = [ cc, tcc + pl ];
    gc = [ gc, tgc + pl ];
    ff = [ff, f_app{i}(so1)]; 
    tt = [tt, f_tra{i}(so1)]; 
    gn = [gn, gapn{i}(so2)]; 
    gt = [gt, gapt{i}(so2)]; 
    gw = [gw, gapw{i}(so2)]; 
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


  if mel_b
    t = tiledlayout(3,1,'TileSpacing','Compact');
  else
    t = tiledlayout(2,1,'TileSpacing','Compact');
  end

  nexttile

  plot(gc,gn,'color',[0.1,0.1,0.6],'linewidth',3,'DisplayName','normal gap')
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
  plot(gc,gt,'color',[0.9,.2,0.3],'linewidth',3,'DisplayName','tangential disp')
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

  if mel_b
    nexttile
    plot(gc,(gw+gn),'color',[0.2,0.1,0.3],'linewidth',3,'DisplayName','rift wall distance')
    hold on
    yl = ylim();
    for i = 1: length(inters)
      plot([inters(i),inters(i)],yl,'c-','linewidth',0.8,'Color',[0.5,0.3,0.7,0.2]);
    end
    ylim(yl);
    xlabel('distance along crack (m)')
    ylabel('wall to wall "distance" (m)')
    tstr = ['"melange" distance between rift walls ',num2str(kk)];
    title(tstr);
  end

  nstr = ['rift',num2str(kk),'_gaps_step',num2str(stepnum)];
  print([results_path,'/',nstr],'-dpng','-r300')

end

