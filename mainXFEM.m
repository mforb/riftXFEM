function [Knum,Theta,xCrk,stop] = mainXFEM(xCrk,npas,delta_inc)

%-- Declare global variables here global elemType stressState typeCrack global L D E nu C P sigmato
global numcrack xCr deltaInc numstep
global plotmesh plotNode plothelp plotiter
global node element numnode numelem bcNodes edgNodes typeProblem elemType
global penalty fixedF contact melange Kpen rift_wall_pressure xM melangeforce stabilize
global epsilon loadstress
global results_path
global fmesh
global output_file
global Hidden zoom_dim 
global wall_int skip_branch skip_vertex modpen modocean
global melange_stab crack_load

output_file = fopen([results_path,'/output.log'],'w')
if ~isfield(xCrk,'tip')
  for i =1:size(xCrk,2)
    xCrk(i).tip = [1,1];
  end
  fprintf(output_file,'No tip field for xCr defined, all tips have been activated\n')
  disp('No tip field for xCr defined, all tips have been activated') ;
end
if isempty(penalty)
  penalty = 0 ;
end
if isempty(contact)
  contact = 0;
end
if isempty(melangeforce)
  melangeforce = 0;
end
if isempty(melange)
  melange = 0;
end
if isempty(wall_int)
  wall_int = 2;
end
if isempty(skip_branch)
  skip_branch = 0;
end
if isempty(skip_vertex)
  skip_vertex = 0;
end
if isempty(modpen)
  modpen = 0;
end
if isempty(modocean)
  modpen = 0;
end
if isempty(melange_stab)
  melange_stab = 0;
end
nitsche = 0;
plot_stresses = 0;
TR = triangulation(element,node);

Knum = [ ] ; Theta = [ ] ;
enrDomain = [ ] ; tipElem = [ ] ; splitElem = [ ] ; vertexElem = [ ] ; cornerElem = [];
%loop over number of steps of crack growth
for ipas = 1:npas
    cgrow = [num2str(toc),'    Crack growth number     ',num2str(ipas)];
    disp(cgrow) ;
    fprintf(output_file,[cgrow,'\n'])
    fprintf(output_file,'---------------------------------------------------\n')
    fprintf(output_file,'---------------------------------------------------\n')
    nstr = ['Number of elements: ',num2str(numelem)];
    fprintf(output_file,[nstr,'\n'])
    

    disp([num2str(toc),'    Crack processing']) ;
    %find elements within a small region
    [enrDomain] = crackDetect(xCrk,ipas,tipElem,splitElem,vertexElem,cornerElem,enrDomain) ;

    %find type of element: tip, split, vertex
    [typeElem,elemCrk,tipElem,splitElem,vertexElem,cornerElem,tangentElem,xTip,xVertex,enrichNode,crackNode,xCrk] = nnodeDetect(xCrk,enrDomain) ;
    % the crack can be slightly modified in cases where crack crosses an element twice
    % if there are any tangent elements
    if ~isempty(tangentElem)
      [nodeTanfix] = f_tangent_iso_node(tangentElem,splitElem,crackNode);
      tan_info = [' CRACK NODES :  ',num2str(length(crackNode)),' crack nodes, ',num2str(length(tangentElem)),' tangent elements, requiring ', num2str(length(nodeTanfix)),' fixed nodes\n'];
      fprintf(output_file,tan_info)
    else
      %tan_element = [];
      %tan_elemCrk = [];
      nodeTanfix = [];
    end
      
    
      %warning('Phantom nodes introduced to account for crack going through nodes')
      %[enrichNode, n_red] =  f_phantomNode(crackNode,elemCrk,splitElem,tipElem,vertexElem,enrichNode) ;
    %else
      %n_red = 0;
    %end 


    



        %plot enriched nodes
    figure(fmesh);
    hold on
    if( strcmp(plotmesh,'YES') )
        for k=1:size(xCr,2)
            for kj = 1:size(xCr(k).coor,1)-1
                cr = plot(xCr(k).coor(kj:kj+1,1),xCr(k).coor(kj:kj+1,2),'r-') ;
                set(cr,'LineWidth',3);
            end
            for kj = 1:size(xCr(k).coor,1)
                plot(xCr(k).coor(kj,1),xCr(k).coor(kj,2),'ro',...
                    'MarkerFaceColor',[.49 1 .63],'MarkerSize',5);
            end
            %corner_on_crack
            split_nodes = find(enrichNode(:,k) == 2);
            tip_nodes   = find(enrichNode(:,k) == 1);
            n1 = plot(node(split_nodes,1),node(split_nodes,2),'r*');
            n2 = plot(node(tip_nodes,1),node(tip_nodes,2),'rs');
            set(n1,'MarkerSize',15);
            set(n2,'MarkerSize',15);
        end
    end
    print([results_path,'/mesh_crack_',num2str(ipas)],'-dpng','-r500')

    
   
    %initialize stiffness matrix, force vector
    %each split node is enriched by ONE function, H(x)
    %each tip node is enriched by FOUR functions, B_i(x), i=1,2,3,4
    %total dof = numnode*ndof + numsplitnode*1*ndof + numtipnode*4*ndof
    split = 0 ; tip = 0 ;
    for k=1:size(xCrk,2)
        split = split + size(find(enrichNode(:,k) == 2), 1) ;
        tip = tip + size(find(enrichNode(:,k) == 1),1 ) ;
    end

    ndof = 2 ;
    %totalUnknown = numnode*ndof + split*1*ndof + tip*4*ndof - n_red*1*ndof ;
    totalUnknown = numnode*ndof + split*1*ndof + tip*4*ndof ;

    K = sparse(totalUnknown,totalUnknown) ;
    F = sparse(totalUnknown,1) ;


    %stiffness matrix computation
    pos = posi(xCrk,numnode,enrichNode,crackNode) ;

    disp([num2str(toc),'    Stiffness Matrix Computation']) ;

    if ~exist('loadstress') || ~strcmp(loadstress,'y')

    [K] = KmatXFEM(enrichNode,elemCrk,typeElem,xTip,xVertex,...
        splitElem,tipElem,vertexElem,cornerElem,crackNode,pos,xCrk,K) ;

    else

      [K,F] = KmatXFEM3(enrichNode,elemCrk,typeElem,xTip,xVertex,...
        splitElem,tipElem,vertexElem,cornerElem,crackNode,enrDomain,pos,xCrk,K,F) ;
    end
    [F] = ForceVector(F) ;

    
    if exist('rift_wall_pressure') & rift_wall_pressure
      disp([num2str(toc),'    Applying ocean pressure imbalance']) ;
      Ft = F;
      [F,elemForce] = f_apply_ocean_pressure(enrichNode,elemCrk,typeElem,xTip,xVertex,...
        splitElem,tipElem,vertexElem,cornerElem,crackNode,enrDomain,[],pos,xCrk,F) ;
    elseif ~isempty(crack_load) & crack_load ~= 0
      disp([num2str(toc),'    Applying load on crack lips']) ;
      [F,elemForce] = f_apply_crack_load(enrichNode,elemCrk,typeElem,xTip,xVertex,...
        splitElem,tipElem,vertexElem,cornerElem,crackNode,enrDomain,[],pos,xCrk,F,crack_load) ;
    else
      elemForce = zeros(2,size(element,1),wall_int*2); % 2 potential segments, all elements, int points * 2 for normal and tangential
    end
    %keyboard

    if exist('melange') & melange 
      if ~isfield(xCrk,'melange') 
        np = size(xM.coor,1);
        if np < 4
          warning('Not enough coordinates in the crack for generic melange aproach')
        else 
          xM.melange = ones(1,np-1);
          xM.melange(1) = 0;
          xM.melange(end) = 0;
          xM.width = ones(100,np-1);
        end
      end
      Kt = sparse(zeros(size(K)));
      [Kt,nodeTanfix] = KmatMELAN(enrichNode,elemCrk,typeElem,xVertex,xTip,...
        splitElem,tipElem,vertexElem,cornerElem,tangentElem,crackNode,pos,xM,xCrk,Kt,nodeTanfix) ;
      K = K + Kt;
      %keyboard
    end

    if ~isempty(nodeTanfix)
      % these are dof that are only connected to crack in a tangent element (therefore enrichment is unconstrained). 
      inds = [];
      figure(1)
      for i = 1:length(nodeTanfix)
        in = [ 2*pos(nodeTanfix(i))-1, 2*pos(nodeTanfix(i))]; 
        inds = [inds, in];
        warning(['fixing tangent dof : ', num2str(in(1)), ' and ', num2str(in(2))] )
        plot(node(nodeTanfix(i),1),node(nodeTanfix(i),2),'sk','linestyle','none','markersize',4)
        
      end
      K(inds,:) = 0;
      K(:,inds) = 0;
      F(inds) = 0;
      for n = 1:length(inds)
        K(inds(n),inds(n)) = 1;
      end
    end



    kk1 = K ;
    %----- Imposing Essential boundary conditions
    disp([num2str(toc),'    Imposing Essential boundary conditions'])
    
    if (strcmp(typeProblem,'eCrkTen') || strcmp(typeProblem,'dispFriction'))
        
        dispNodes = unique([bcNodes{1}]) ;
        bcdof = [ ]; bcval = [ ];
        for i=1:length(dispNodes)
            bcdof = [bcdof 2*dispNodes(i)] ;
            bcval = [bcval 0] ;
        end
        bcdof = [bcdof 2*dispNodes(1)-1] ;
        bcval = [bcval 0];
        if strcmp(typeProblem,'dispFriction')
          dispNodes = unique([bcNodes{3}]) ; % the top nodes
          for i=1:length(dispNodes)
              bcdof = [bcdof 2*dispNodes(i)] ;
              bcval = [bcval -0.1] ;
          end
        end
    elseif strcmp(typeProblem,'eCrkTen2')
        dispNodes = unique([bcNodes{1}]) ;
        bcdof = [ ]; bcval = [ ];
        for i=1:length(dispNodes)
            bcdof = [bcdof 2*dispNodes(i)-1 2*dispNodes(i)] ;
            bcval = [bcval 0 0] ;
        end
    elseif strcmp(typeProblem,'centre')
        dispNodes = unique([bcNodes{3}]) ;
        bcdof = [ ]; bcval = [ ];
        for i=1:length(dispNodes)
            bcdof = [bcdof 2*dispNodes(i)] ;
            bcval = [bcval 0] ;
        end
        dispNodes = unique([bcNodes{4}]) ;
        for i=1:length(dispNodes)
            bcdof = [bcdof 2*dispNodes(i)-1] ;
            bcval = [bcval 0 ] ;
        end
    elseif strcmp(typeProblem,'ISSM')
        dispNodes = [bcNodes{4}];
        bcdof = [ ]; bcval = [ ];
        for i=1:length(dispNodes)
            bcdof = [bcdof 2*dispNodes(i)-1 2*dispNodes(i)] ;
            bcval = [bcval 0 0] ;
        end
    elseif strcmp(typeProblem,'Test')
      dispNodes = unique([bcNodes{4}]);
      bcdof = [2*dispNodes(end)-1 2*dispNodes(end)]; 
      bcval = [ 0 0 ];
    else
        %boundary condition
        cracklength = 100 ;
        
        xcTip = xCr(1).coor(2,:) ;
        seg = xCr(1).coor(2,:) - xCr(1).coor(1,:) ;
        dispNodes  = [bcNodes{1};bcNodes{2};bcNodes{3};bcNodes{4}] ;
        dispNodes = unique(dispNodes) ;
        bcdof = [ ]; bcval = [ ] ;
        for i = 1:length(dispNodes)
            bcdof = [bcdof 2*dispNodes(i)-1 2*dispNodes(i)] ;
            x = node(dispNodes(i),:) ;
            [ux,uy] = exact_Griffith(x,E,nu,stressState,P,xcTip,seg,cracklength) ;
            uFixed(i) = ux ;
            vFixed(i) = uy ;
            bcval = [bcval uFixed(i) vFixed(i)] ;
        end
    end
    
    %number of dofs
    numdof = totalUnknown - size(bcdof,2) ;
    
    %apply boundary conditions and solve
    
    [K,F] = feaplyc2(K,F,bcdof,bcval) ;

    %[L,U] = lu(K) ;
    %y = L\F;
    %u = U\y;
    u = K\F;

    fu = full(u);
    Stdux = fu(1:2:2*numnode) ;
    Stduy = fu(2:2:2*numnode) ;
    %[crackLips,flagP] = f_cracklips( u, xCrk, enrDomain, typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);

    [crackLips,flagP,elemGap] = f_find_cracklips( u, xCrk, 1, enrDomain, typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
    if Hidden 
      f = figure('visible','off');
    else
      f = figure();
    end
    hold on
    dfac = 1 ;
    triplot(TR);
    hold on
    f_plotCrack(crackLips,2000,'r-','k-','c--')
    print([results_path,'/crack_walls_before',num2str(ipas)],'-dpng','-r500')
    if ~isempty(zoom_dim)
      xlim(zoom_dim(1,:));
      ylim(zoom_dim(2,:));
      figure_name = ['crack_walls_before_zoom',num2str(ipas)];
      print([results_path,'/',figure_name],'-dpng','-r500')
      %keyboard
    end
    clf();
    trisurf(element,node(:,1),node(:,2),Stduy)
    axis equal; view(2); shading interp; colorbar
    title('Y displacement before Newton solver')
    print([results_path,'/original_ydisp',num2str(ipas)],'-dpng')
    clf();


    if contact & ~flagP
      % first we need to find out if there is any interpenetration
      penalty = 0
      disp([num2str(toc),'    No contact therefore penalty method was not applied'])
    elseif ( contact & flagP ) | melange_stab
      if ~isempty(stabilize) & stabilize
        K_orig = K;
        %Rc1 = rcond(full(K));
        [K] = KmatSTAB(Kpen,enrichNode,crackNode,elemCrk,typeElem,xTip,xVertex,splitElem,tipElem,vertexElem,cornerElem,tangentElem,pos,xCrk,K,u,melange_stab);
        %Rc2 = rcond(full(K));
        %stab_str = ['STABALIZATION: enriched dofs. Rcond is originally ',num2str(Rc1),' and then becomes ',num2str(Rc2),'\n'];
        %fprintf(output_file,stab_str)
        disp([num2str(toc),'    Recalculating u with stabalized K'])
        u = K\F;
        fu2 = full(u);
        Stdux2 = fu2(1:2:2*numnode) ;
        Stduy2 = fu2(2:2:2*numnode) ;
        if Hidden 
          f = figure('visible','off');
        else
          f = figure();
        end
        trisurf(element,node(:,1),node(:,2),Stduy2-Stduy)
        axis equal; view(2); shading interp; colorbar
        title('Y displacement change from stabilizing')
        print([results_path,'/stab_ydiff_',num2str(ipas)],'-dpng')
        if ~isempty(zoom_dim)
          xlim(zoom_dim(1,:));
          ylim(zoom_dim(2,:));
          figure_name = ['stab_ydiff_zoom_',num2str(ipas)];
          print([results_path,'/',figure_name],'-dpng','-r500')
        end
        figure(f)
        clf()
        [crackLips,flagP] = f_find_cracklips( u, xCrk, 1, enrDomain, typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
        dfac = 1 ;
        plotMesh(node+dfac*[Stdux, Stduy],element,elemType,'b-',plotNode,f)
        f_plotCrack(crackLips,20,'r-','k-','c--')
        if ~isempty(zoom_dim)
          xlim(zoom_dim(1,:));
          ylim(zoom_dim(2,:));
        end
        print([results_path,'/crack_walls_stab',num2str(ipas)],'-dpng','-r500')
        figure(f)
        clf()
        trisurf(element,node(:,1),node(:,2),Stduy2)
        axis equal; view(2); shading interp; colorbar
        title('Y displacement after stabilizing')
        print([results_path,'/stab_ydisp_',num2str(ipas)],'-dpng')
        if ~isempty(zoom_dim)
          xlim(zoom_dim(1,:));
          ylim(zoom_dim(2,:));
          figure_name = ['stab_ydisp_zoom_',num2str(ipas)];
          print([results_path,'/',figure_name],'-dpng','-r300')
        end
        clf()
      end
      penalty = 1;
    end

    if melangeforce
      penalty = 1;
    end


    if penalty 
      elemForce_orig = elemForce;
      elemForce = zeros(size(elemForce));
      tol1 = 1e-10;
      %tol2 = 1e-8;
      tol2 = 1e-2;
      cont = 1
      Du = zeros(size(u));
      Fext = F
      %if melangeforce
        %u = u/2;
      %end
      nu = 1
      while 1 
        disp(['Newton step ',num2str(cont)])
        disp(['---------------------------------------------'])
        Fint = K*u;
        if nitsche
          [KT,Gint,elemForce] = KmatNITSCHE(Kpen,enrichNode,crackNode,elemCrk,typeElem,xTip,xVertex,splitElem,tipElem,vertexElem,cornerElem,tangentElem,elemForce,pos,xCrk,xM,K,u);
        %[KT] = KmatSTAB(Kpen,enrichNode,crackNode,elemCrk,typeElem,xTip,xVertex,splitElem,tipElem,vertexElem,cornerElem,tangentElem,pos,xCrk,KT,u);
        else 
          [KT,Gint,elemForce,ratc] = KTmatXFEM(Kpen,enrichNode,crackNode,elemCrk,typeElem,xTip,xVertex,splitElem,tipElem,vertexElem,cornerElem,tangentElem,elemForce,pos,xCrk,xM,K,u);
        end
        Res  = Fext - Fint - Gint ;
        nr = norm(Res,2);
        if cont == 1
          nr0  = nr;
        end
        rnr = nr/nr0;
        disp(['L2 norm of the residual, R =  ',num2str(nr)])
        disp(['Relative to R0 : ',num2str(rnr)])
        out_str = ['L2 norm of residual is ',num2str(nr),' ; relative to R0 ',num2str(rnr),'. The contact ratio is: ',num2str(ratc),'\n'];
        fprintf(output_file,out_str)
        
        if rnr < tol1 | nr < tol2
           disp(['Converged at step : ',num2str(cont)])
           break
        %elseif cont > 200
        elseif cont > 100 
           warning(['After, ',num2str(cont),' iterations ||R||/||R0|| is still: ',num2str(rnr)])
           break
        end
        u = u + KT\Res;
        cont = cont + 1;

        if plothelp | plotiter 
          fu = full(u);
          Stdux = fu(1:2:2*numnode) ;
          Stduy = fu(2:2:2*numnode) ;
          [crackLips,flagP] = f_find_cracklips( u, xCrk, 1, enrDomain, typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
          f = figure('visible','on');
          clf
          hold on
          dfac = 1 ;
          plotMesh(node+dfac*[Stdux, Stduy],element,elemType,'b-',plotNode,f)
          f_plotCrack(crackLips,1,'r-','k-','m--')
          print(['crack_iter',num2str(cont)],'-dpng','-r500')
        end
      end
      [crackLips,flagP,elemGap] = f_find_cracklips( u, xCrk, 1, enrDomain, typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
      fu = full(u);
      Stdux = fu(1:2:2*numnode) ;
      Stduy = fu(2:2:2*numnode) ;
      try
        [inters,gn_inters] = f_plot_wall_forces(u,xCrk,enrDomain,typeElem,elemForce,elemGap,elemCrk,splitElem,vertexElem,tipElem,ipas);
        elemForce = elemForce + elemForce_orig;
        f_plot_wall_forces(u,xCrk,enrDomain,typeElem,elemForce,elemGap,elemCrk,splitElem,vertexElem,tipElem,ipas+100);
        %f_plot_wall_forces(u,xCrk,enrDomain,typeElem,elemForce_orig,elemGap,elemCrk,splitElem,vertexElem,tipElem,ipas+200)
      catch
        keyboard
      end

%     
  %     % plot displacement contour
      if Hidden 
        f = figure('visible','off');
      else
        f = figure();
      end
      trisurf(element,node(:,1),node(:,2),Stduy)
      axis equal; view(2); shading interp; colorbar
      title(['Y displacement after Newton solver'])
      print([results_path,'/after_ydisp'],'-dpng')
      clf(f)
       
       %save('test.mat','K','F','u')

      if Hidden 
        f = figure('visible','off');
      else
        f = figure();
      end
      hold on
      dfac = 1 ;
      %plotMesh(node+dfac*[Stdux, Stduy],element,elemType,'b-',plotNode,f)
      triplot(TR);
      f_plotCrack(crackLips,200,'r-','k-','c--')
      %keyboard
      print([results_path,'/crack_walls_after',num2str(ipas)],'-dpng','-r500')
      if ~isempty(zoom_dim)
        xlim(zoom_dim(1,:));
        ylim(zoom_dim(2,:));
        figure_name = ['crack_walls_after_zoom',num2str(ipas)];
        print([results_path,'/',figure_name],'-dpng','-r500')
        %keyboard
      end
      clf(f)
    end
    if melange
      [inters,gn_inters] = f_plot_wall_forces(u,xCrk,enrDomain,typeElem,elemForce,elemGap,elemCrk,splitElem,vertexElem,tipElem,ipas);
    end

    if plot_stresses 
      if strcmp(elemType,'Q4')
         plotFieldXfem(xCrk,pos,enrichNode,u,...
             elemCrk,vertexELem,splitElem,tipElem,xVertex,xTip,typeElem,ipas);    
      else
        plotFieldXfemT3(xCrk,pos,enrichNode,crackNode,u,...
          elemCrk,vertexElem,cornerElem,splitElem,tipElem,xVertex,xTip,typeElem,ipas) ;
        %keyboard
      end
    end

    if Hidden 
      f = figure('visible','off');
    else
      f = figure();
    end
    hold on
    dfac = 4000 ;
    plotMesh(node+dfac*[Stdux, Stduy],element,elemType,'b-',plotNode,f)
     %plotMesh(node+dfac*[uxAna uyAna],element,elemType,'r-',plotNode)
    figure_name = ['Disp_fact_',num2str(dfac),'_',num2str(ipas)];
    print([results_path,'/',figure_name],'-dpng','-r300')
    clf(f)

    if ~exist('gn_inters')
      gn_inters = zeros(length(xCrk.coor),1);
    end
   
    [Knum,Theta,xCrk,stop] = KcalJint(xCrk,...
        typeElem,enrDomain,elemCrk,enrichNode,crackNode,xVertex,...
        vertexElem,pos,u,ipas,delta_inc,Knum,Theta,...
        tipElem,splitElem,cornerElem,elemForce,gn_inters) ;

    %keyboard
    var_name = [results_path,'/crack',num2str(ipas),'.mat'];
    save(var_name,'xCrk','Knum','Theta','u','element','node','pos','enrichNode','crackNode','elemCrk','vertexElem','cornerElem','splitElem','tipElem','xVertex','xTip','typeElem','bcNodes','elemForce','elemGap');
    if stop
      break;
    end
end
fclose(output_file);
