function [inters,ylg,yls] = f_plot_wf( u,xCrk,enrDomain,typeElem,elem_force,elem_gap,elem_crk,split_elem,vertex_elem,xVertex,tip_elem,xTip,crack_node,type_elem,enrich_node,pos,stepnum,varargin)
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
intType = 'TRIANGULAR';
corner = [1 2 3 1] ;
nnode = [0 0;1 0;0 1] ;

if isempty(wall_int)
  wall_int = 2;
end

if isempty(melange)  | melange == 0 
  mel_b = 0;
else
  mel_b = 1
end

if nargin == 19
  ylg = varargin{1}; 
  yls = varargin{2}; 
else
  ylg = []; 
  yls = []; 
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
  s_n   = {[]};
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
        s_n{i}  = [];
        gapn{i} = [];
        gapt{i} = [];
        gapw{i}  = [];
      end

      xseg = [xCrk(kk).coor(i,:),xCrk(kk).coor(i+1,:)];
      [lseg,~,~,~,~,~] = f_segment_dist(xseg);
      [flag1,flag2,~] = crack_interact_element(xseg,iel,[]);
      if flag1
        % local stuff
        seg = xCrk(kk).coor(i,:) - xCrk(kk).coor(i+1,:); 
        alpha = atan2(seg(2),seg(1));
        QT  =  [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];       
        [ap,apg] = f_crack_wall(iel,nnode,corner,tip_elem,vertex_elem,elem_crk,xTip,xVertex,crack_node);
        ap = f_align_lp_gc(ap,[apg(1,:),apg(end,:)],sctr);
        seg = 1;
        p = ap(seg:seg+1,:);
        sv = [];
        for k_in = 1:length(Q)
          [Np,dNdxp]=lagrange_basis('L2',Q(k_in));
          gpt = Np'*p;
          s = f_calc_stress(gpt,iel,u,type_elem,enrich_node,elem_crk,xVertex,xTip,crack_node,pos,kk,QT);
          sv = [sv, s(2,2)];
        end
        % global stuff
        [l,~,~,~,~,~] = f_segment_dist(segment);
        ccoord = ((l/2).*(Q+1))';
        [ccoord,so] = sort(ccoord);
        ccoord = [0,ccoord]; % for the ends
        gcoord = [0,l];
        ind = 4;
        fp = elem_force(1,iel,1:2:wall_int*2); fp=fp(so); fp = reshape(fp,1,wall_int);
        ft = elem_force(1,iel,2:2:wall_int*2); ft = ft(so); ft = reshape(ft,1,wall_int);
        sv = sv(so); sv = reshape(sv,1,wall_int);

        if ~points_same_2d(segment(1:2),xseg(1:2))
          ccoord = lseg - l + ccoord; 
          ccoord(1) = lseg; % the first coordinate is still the end
          gcoord = lseg - gcoord;
          ind = 2;
        end
        crack_coord{i} = [crack_coord{i}, ccoord];
        f_app{i} = [ f_app{i},0,fp ];
        f_tra{i} = [ f_tra{i},0,ft ];
        s_n{i} = [ s_n{i},0,ft ];
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
        % local stuff
        seg = xCrk(kk).coor(i,:) - xCrk(kk).coor(i+1,:); 
        alpha = atan2(seg(2),seg(1));
        QT  =  [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];       
        [ap,apg] = f_crack_wall(iel,nnode,corner,tip_elem,vertex_elem,elem_crk,xTip,xVertex,crack_node);
        ap = f_align_lp_gc(ap,[apg(1,:),apg(end,:)],sctr);
        seg = 1;
        p = ap(seg:seg+1,:);
        sv = [];
        for k_in = 1:length(Q)
          [Np,dNdxp]=lagrange_basis('L2',Q(k_in));
          gpt = Np'*p;
          s = f_calc_stress(gpt,iel,u,type_elem,enrich_node,elem_crk,xVertex,xTip,crack_node,pos,kk,QT);
          sv = [sv, s(2,2)];
        end
        % global stuff
        nseg = [xseg(1:2),segment(1:2)];
        [nl,~,~,~,~,~] = f_segment_dist(nseg);
        ccoord = (l/2).*(Q+1)'+nl;
        ccoord = ccoord(so);
        gc     = [nl,nl+l];
        fp = elem_force(1,iel,1:2:wall_int*2); fp = fp(so); fp = reshape(fp,1,wall_int);
        ft = elem_force(1,iel,2:2:wall_int*2); ft = ft(so); ft = reshape(ft,1,wall_int);
        sv = sv(so); sv = reshape(sv,1,wall_int);
        crack_coord{i} = [crack_coord{i}, ccoord];
        crack_gapc{i}   = [crack_gapc{i},gc];
        f_app{i} =  [f_app{i},fp];
        f_tra{i} = [f_tra{i},ft];
        s_n{i} = [s_n{i},sv];
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
        % local stuff
        seg = xCrk(kk).coor(i,:) - xCrk(kk).coor(i+1,:); 
        alpha = atan2(seg(2),seg(1));
        QT  =  [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];       
        [ap,apg] = f_crack_wall(iel,nnode,corner,tip_elem,vertex_elem,elem_crk,xTip,xVertex,crack_node);
        ap = f_align_lp_gc(ap,[apg(1,:),apg(end,:)],sctr);
        seg = 1;
        p = ap(seg:seg+1,:);
        seg = 2;
        p2 = ap(seg:seg+1,:);
        sv1 = [];
        sv2 = [];
        for k_in = 1:length(Q)
          [Np,dNdxp]=lagrange_basis('L2',Q(k_in));
          gpt = Np'*p;
          gpt2 = Np'*p2;
          s = f_calc_stress(gpt,iel,u,type_elem,enrich_node,elem_crk,xVertex,xTip,crack_node,pos,kk,QT);
          sv1 = [sv1, s(2,2)];
          s = f_calc_stress(gpt2,iel,u,type_elem,enrich_node,elem_crk,xVertex,xTip,crack_node,pos,kk,QT);
          sv2 = [sv2, s(2,2)];
        end

        % global stuff
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
        sv1 = sv1(so); sv1 = reshape(sv1,1,wall_int);
        %fp1 = (fp1)/d1;
        %ft1 = (ft1)/d1;
        fp2 = elem_force(2,iel,1:2:wall_int*2); fp2 = fp2(so); fp2 = reshape(fp2,1,wall_int);
        ft2 = elem_force(2,iel,2:2:wall_int*2); ft2 = ft2(so); ft2 = reshape(ft2,1,wall_int);
        sv2 = sv2(so); sv2 = reshape(sv2,1,wall_int);
        %fp2 = (fp2)/d1;
        %ft2 = (ft2)/d2;
        crack_coord{i} = [ crack_coord{i}, ccoord1 ];
        f_app{i} = [ f_app{i}, fp1 ];
        f_tra{i} = [ f_tra{i}, ft1 ];
        s_n{i} = [ s_n{i}, sv1 ];
        crack_gapc{i} = [ crack_gapc{i}, gc1 ];
        gapn{i} = [ gapn{i}, elem_gap(iel,[2,4]) ];
        gapt{i} = [ gapt{i}, elem_gap(iel,[1,3]) ];
        gapw{i} = [ gapw{i}, width,width ];
        crack_coord{i+1} = [ ccoord2, crack_coord{i+1} ];
        f_app{i+1} = [ fp2, f_app{i+1}];
        f_tra{i+1} = [ ft2, f_tra{i+1}];
        s_n{i+1} = [ sv2, s_n{i+1}];
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
  ss = [];
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
    ss = [ss, s_n{i}(so1)]; 
    gn = [gn, gapn{i}(so2)]; 
    gt = [gt, gapt{i}(so2)]; 
    gw = [gw, gapw{i}(so2)]; 
    pl = pl + lseg;
    inters = [inters, pl];
  end
  xl = [ 0 inters(end)/1000];
  inters(end) = [];

  t = tiledlayout(2,1,'TileSpacing','Compact');
  nexttile

  plot(cc/1000,ff,'linewidth',3,'DisplayName','normal sym force')
  yl = ylim();
  hold on
  for i = 1: length(inters)
    pi(i) = plot([inters(i),inters(i)]/1000,yl,'c-','linewidth',2,'Color',[0.6,0.3,0.5,0.2]);
  end
  ylim(yl);
  xlim(xl);
  xlabel('distance along crack')
  ylabel('force')
  tstr = ['normal forces acting on rift ',num2str(kk)];
  title(tstr);

  nexttile
  plot(cc/1000,tt,'color',[1,0.1,0],'linewidth',2,'DisplayName','tangential force')
  hold on
  yl = ylim();
  for i = 1: length(inters)
    pi(i) = plot([inters(i),inters(i)]/1000,yl,'c-','linewidth',2,'Color',[0.6,0.3,0.5,0.2]);
  end
  ylim(yl);
  xlim(xl);
  xlabel('distance along crack')
  ylabel('force')
  tstr = ['tangential forces acting on rift ',num2str(kk)];
  title(tstr);
  nstr = ['rift',num2str(kk),'forces_step',num2str(stepnum)];
  print([results_path,'/',nstr],'-dpng','-r300')
  clf

  if Hidden
    f = figure('visible','off');
  else
    f = figure();
  end
  plot(cc/1000,ss,'color',[0.04,0.04,0.04],'linewidth',3,'DisplayName','normal sym force')
  yl = ylim();
  hold on
  for i = 1: length(inters)
    pi(i) = plot([inters(i),inters(i)]/1000,yl,'c-','linewidth',2,'Color',[0.6,0.3,0.5,0.2]);
  end
  ylim(yl);
  xlim(xl);
  xlabel('distance along crack')
  ylabel('crack normal stress')
  nstr = ['rift',num2str(kk),'_Nstress_step',num2str(stepnum)];
  print([results_path,'/',nstr],'-dpng','-r300')
  
  if Hidden
    f = figure('visible','off');
  else
    f = figure();
  end
  pg = plot(gc/1000,gn,'color',[0.04,0.04,0.04],'linewidth',2,'DisplayName','normal gap')
  yl = ylim();
  hold on
  if isempty(ylg)
    ylg = yl;
  else
    yl = ylg;
  end
  for i = 1: length(inters)
    pi(i) = plot([inters(i),inters(i)]/1000,yl,'c-','linewidth',2,'Color',[0.6,0.3,0.5,0.2]);
  end
  ylim(yl);
  xlim(xl);
  %ax.XAxisLocation = 'origin';
  grid on
  xlabel('Distance along rift (km)')
  ylabel('Gap (m)')
  tstr = ['normal gap along rift ',num2str(kk)];
  %title(tstr);
  nstr = ['rift',num2str(kk),'_gap_',num2str(stepnum)];
  print([results_path,'/',nstr],'-dpng','-r300')
  if stepnum == 1
    delete(pi)
    pi = []
    yl = [-.55,.55];
    for i = 1: length(inters)
      pi(i) = plot([inters(i),inters(i)]/1000,yl,'c-','linewidth',2,'Color',[0.6,0.3,0.5,0.2]);
    end
    ylim(yl)
    set(pg,'linewidth',3);
    nstr = ['rift',num2str(kk),'_gap_normalized'];
    print([results_path,'/',nstr],'-dpng','-r300')
  end

  clf

  plot(gc/1000,gt,'color',[0.04,0.04,0.04],'linewidth',3,'DisplayName','tangential disp')
  hold on
  yl = ylim();
  if isempty(yls)
    yls = yl;
  else
    yl = yls;
  end
  for i = 1: length(inters)
    pi(i) = plot([inters(i),inters(i)]/1000,yl,'c-','linewidth',2,'Color',[0.6,0.3,0.5,0.2]);
  end
  xlim(xl);
  ylim(yl);
  xlabel('distance along crack (m)')
  ylabel('slip (m)')
  tstr = ['tangential displacement along rift ',num2str(kk)];
  nstr = ['rift',num2str(kk),'_slip_',num2str(stepnum)];
  print([results_path,'/',nstr],'-dpng','-r300')
  clf

  if mel_b
    plot(gc/2,(gw+gn),'color',[0.04,0.04,0.04],'linewidth',2,'DisplayName','rift wall distance')
    hold on
    yl = ylim();
    for i = 1: length(inters)
      pi(i) = plot([inters(i),inters(i)]/1000,yl,'c-','linewidth',2,'Color',[0.6,0.3,0.5,0.2]);
    end
    ylim(yl);
    xlabel('distance along crack (m)')
    ylabel('wall to wall "distance" (m)')
    tstr = ['"melange" distance between rift walls ',num2str(kk)];
    nstr = ['rift',num2str(kk),'_walld_',num2str(stepnum)];
    print([results_path,'/',nstr],'-dpng','-r300')
  end

end

