function [Knum,Theta,xCrk] = mainXFEM(xCrk,npas,delta_inc)

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
global wall_int

output_file = fopen([results_path,'/output.log'],'w')
if ~isfield(xCr,'tip')
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
if isempty(wall_int)
  wall_int = 2;
end

Knum = [ ] ; Theta = [ ] ;
enrDomain = [ ] ; tipElem = [ ] ; splitElem = [ ] ; vertexElem = [ ] ; cornerElem = [];
%loop over number of steps of crack growth
for ipas = 1:npas
    cgrow = [num2str(toc),'    Crack growth number     ',num2str(ipas)];
    disp(cgrow) ;
    fprintf(output_file,[cgrow,'\n'])
    fprintf(output_file,'---------------------------------------------------\n')
    fprintf(output_file,'---------------------------------------------------\n')
    

    disp([num2str(toc),'    Crack processing']) ;
    %find elements within a small region
    [enrDomain] = crackDetect(xCrk,ipas,tipElem,splitElem,vertexElem,cornerElem,enrDomain) ;

    %find type of element: tip, split, vertex
    [typeElem,elemCrk,tipElem,splitElem,vertexElem,cornerElem,tangentElem,xTip,xVertex,enrichNode,crackNode] = nnodeDetect(xCrk,enrDomain) ;
    % if there are any tangent elements
    if ~isempty(tangentElem)
      [nodeTanfix] = f_tangent_iso_node(tangentElem,crackNode);
      tan_info = [' CRACK NODES :  ',num2str(length(crackNode)),' crack nodes, ',num2str(length(tangentElem)),' tangent elements, requiring ', num2str(length(nodeTanfix)),' fixed nodes\n'];
      fprintf(output_file,cgrow)
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
    else
      elemForce = zeros(2,size(element,1),wall_int*2); % 2 potential segments, all elements, int points * 2 for normal and tangential
    end

    if exist('melange') & melange 
      if ~exist('xM') 
        xM = xCrk; % this means that all of the crack including the tip elements will be infilled with melange
        np = size(xM.coor,1);
        if np < 4
          warning('Not enough coordinates in the crack for generic melange aproach')
        else 
          xM.melange = ones(1,np-1);
          xM.melange(1) = 0;
          xM.melange(end) = 0;
          xM.width = ones(500,np-1);
        end
      end
      Kt = K;
      Kt = K - Kt;
      [Kt] = KmatMELAN(enrichNode,elemCrk,typeElem,xVertex,xTip,...
        splitElem,tipElem,vertexElem,cornerElem,tangentElem,crackNode,pos,xM,xCrk,Kt) ;
      K = K + Kt;
      %keyboard
    end

    if ~isempty(nodeTanfix)
      % these are dof that are only connected to crack in a tangent element (therefore enrichment is unconstrained). 
      inds = [];
      for i = 1:length(nodeTanfix)
        in = [ 2*pos(nodeTanfix(i))-1, 2*pos(nodeTanfix(i))]; 
        inds = [inds, in];
        warning(['fixing tangent dof : ', num2str(in(1)), ' and ', num2str(in(2))] )
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

    %if exist('fixedF') & ~isempty(fixedF)
      %[crackLips,flagP] = f_cracklips( zeros(totalUnknown,1), xCr, enrDomain, typeElem, elemCrk, xTip, xVertex, enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
      %Fcrack = zeros(size(F));
      %Fcrack = f_crackforce_fixed(Fcrack,fixedF,crackLips,xCr,elemCrk,xTip,pos,typeElem, enrichNode,splitElem,vertexElem,tipElem);
      %keyboard
      %F = F + Fcrack;
    %end

    %[L,U] = lu(K) ;
    %y = L\F;
    %u = U\y;
    u = K\F;

    fu = full(u);
    Stdux = fu(1:2:2*numnode) ;
    Stduy = fu(2:2:2*numnode) ;
    [crackLips,flagP] = f_cracklips( u, xCrk, enrDomain, typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);

    if Hidden 
      f = figure('visible','off');
    else
      f = figure();
    end
    hold on
    dfac = 1 ;
    plotMesh(node+dfac*[Stdux, Stduy],element,elemType,'b-',plotNode,f)
    figure(f);
    hold on
    f_plotCrack2(crackLips,20,'r-','k-','c--')
    print([results_path,'/crack_walls_before',num2str(ipas)],'-dpng','-r300')
    clf(f)
    trisurf(element,node(:,1),node(:,2),Stduy)
    axis equal; view(2); shading interp; colorbar
    title('Y displacement before Newton solver')
    print([results_path,'/original_ydisp',num2str(ipas)],'-dpng')
    clf(f)


    if contact & ~flagP
      % first we need to find out if there is any interpenetration
      contact = 0
      disp([num2str(toc),'    No contact therefore penalty method was not applied'])
    elseif contact & flagP
      if ~isempty(stabilize) & stabilize
        [K] = KmatSTAB(Kpen,enrichNode,crackNode,elemCrk,typeElem,xTip,xVertex,splitElem,tipElem,vertexElem,cornerElem,tangentElem,pos,xCrk,K,u);
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
          print([results_path,'/',figure_name],'-dpng','-r300')
        end
        clf(f)
        [crackLips,flagP] = f_cracklips( u, xCrk, enrDomain, typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
        dfac = 1 ;
        plotMesh(node+dfac*[Stdux, Stduy],element,elemType,'b-',plotNode,f)
        f_plotCrack2(crackLips,20,'r-','k-','c--')
        print([results_path,'/crack_walls_stab',num2str(ipas)],'-dpng','-r300')
        clf(f)
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
        clf(f)
      end
      penalty = 1
    end

    if melangeforce
      penalty = 1
    end


    if penalty 
      elemForce_orig = elemForce;
      elemForce = zeros(size(elemForce));
      tol1 = 1e-15;
      tol2 = 1e-13;
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
        [KT,Gint,elemForce] = KTmatXFEM(Kpen,enrichNode,crackNode,elemCrk,typeElem,xTip,xVertex,splitElem,tipElem,vertexElem,cornerElem,tangentElem,elemForce,pos,xCrk,xM,K,u);
        Res  = Fext - Fint - Gint ;
        nr = norm(Res,2);
        if cont == 1
          nr0  = nr;
        end
        rnr = nr/nr0;
        disp(['L2 norm of the residual, R =  ',num2str(nr)])
        disp(['Relative to R0 : ',num2str(rnr)])
        if rnr < tol1 | nr < tol2
           disp(['Converged at step : ',num2str(cont)])
           break
        %elseif cont > 200
        elseif cont > 50
           warning(['After, ',num2str(cont),' iterations ||R||/||R0|| is still: ',num2str(rnr)])
           break
        end
        u = u + KT\Res;
        cont = cont + 1;

        if plothelp | plotiter 
          fu = full(u);
          Stdux = fu(1:2:2*numnode) ;
          Stduy = fu(2:2:2*numnode) ;
          [crackLips,flagP] = f_cracklips( u, xCrk, enrDomain, typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
          f = figure('visible','on');
          clf
          hold on
          dfac = 1 ;
          plotMesh(node+dfac*[Stdux, Stduy],element,elemType,'b-',plotNode,f)
          f_plotCrack(crackLips,1,'r-','k-','m--')
          print(['crack_iter',num2str(cont)],'-dpng','-r300')
        end
      end
      fu = full(u);
      Stdux = fu(1:2:2*numnode) ;
      Stduy = fu(2:2:2*numnode) ;
      %elemForce = elemForce + elemForce_orig;
      f_plot_wall_forces(u,xCrk,enrDomain,typeElem,elemForce,elemCrk,splitElem,vertexElem,tipElem,ipas)
    end
%     
%     % plot displacement contour
    figure(f);
    trisurf(element,node(:,1),node(:,2),Stduy)
    axis equal; view(2); shading interp; colorbar
    title(['Y displacement after Newton solver'])
    print([results_path,'/after_ydisp'],'-dpng')
    clf(f)
     
     %save('test.mat','K','F','u')

    [crackLips,flagP] = f_cracklips( u, xCr, enrDomain, typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
    figure(f);
    hold on
    dfac = 1 ;
    plotMesh(node+dfac*[Stdux, Stduy],element,elemType,'b-',plotNode,f)
    f_plotCrack2(crackLips,20,'r-','k-','c--')
    print([results_path,'/crack_walls_after',num2str(ipas)],'-dpng','-r300')
    if ~isempty(zoom_dim)
      xlim(zoom_dim(1,:));
      ylim(zoom_dim(2,:));
      figure_name = ['crack_walls_after_zoom',num2str(ipas)];
      print([results_path,'/',figure_name],'-dpng','-r300')
    end
    clf(f)

    
%     res = [Stdux Stduy] ;
%     
%     
%     % plot displacement contour
%     figure
%     clf
%     trisurf(element,node(:,1),node(:,2),uyAna)
%     axis equal; view(2); shading interp; colorbar
%     title('Displacement from Analytical solution')
%     
%     % plot displacement contour
%     figure
%     clf
%     trisurf(element,node(:,1),node(:,2),(Stduy-uyAna))
%     axis equal; view(2); shading interp; colorbar
%     title('Displacement from XFEM-Analytical solution')
%     
%     
%     plotFieldXfem(xCrk,pos,enrichNode,u,...
%         elemCrk,vertexElem,splitElem,tipElem,xVertex,xTip,typeElem) ;
    
     if strcmp(elemType,'Q4')
        plotFieldXfem(xCrk,pos,enrichNode,u,...
            elemCrk,vertexELem,splitElem,tipElem,xVertex,xTip,typeElem,ipas);    
     else
       plotFieldXfemT3(xCrk,pos,enrichNode,crackNode,u,...
         elemCrk,vertexElem,cornerElem,splitElem,tipElem,xVertex,xTip,typeElem,ipas) ;
     end

     figure(f);
     hold on
     dfac = 50 ;
     plotMesh(node+dfac*[Stdux, Stduy],element,elemType,'b-',plotNode,f)
      %plotMesh(node+dfac*[uxAna uyAna],element,elemType,'r-',plotNode)
     figure_name = ['Disp_fact_',num2str(dfac),'_',num2str(ipas)];
     print([results_path,'/',figure_name],'-dpng','-r300')
     clf(f)


   
    [Knum,Theta,xCrk] = KcalJint(xCrk,...
        typeElem,enrDomain,elemCrk,enrichNode,crackNode,xVertex,...
        vertexElem,pos,u,ipas,delta_inc,Knum,Theta,...
        tipElem,splitElem,cornerElem,elemForce) ;

    %keyboard
    var_name = [results_path,'/crack',num2str(ipas),'.mat'];
    save([results_path,'/crack_disp.mat'],'xCrk','Knum','Theta','u','element','node');
end
fclose(output_file);
