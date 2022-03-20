function [Knum,Theta,xCrk] = mainXFEM(xCrk,npas,delta_inc)

%-- Declare global variables here global elemType stressState typeCrack global L D E nu C P sigmato
global numcrack xCr deltaInc numstep
global plotmesh plotNode
global node element numnode numelem bcNodes edgNodes typeProblem elemType
global penalty fixedF contact melange Kpen
global epsilon loadstress
global results_path
global rift_wall_pressure


if ~exist('penalty')
  penalty = 0 ;
end
if ~exist('contact')
  contact = 0;
elseif contact
  penalty = 1 ; % contact is implemented via the penalty method
end
if ~exist('melange')
  melange = 0
elseif melange
  penalty = 1; % material properties within the crack are implemented via penalty method
end


Knum = [ ] ; Theta = [ ] ;
enrDomain = [ ] ; tipElem = [ ] ; splitElem = [ ] ; vertexElem = [ ] ; cornerElem = [];
%loop over number of steps of crack growth
for ipas = 1:npas
    disp([num2str(toc),'    Crack growth number     ',num2str(ipas)]) ;

    disp([num2str(toc),'    Crack processing']) ;
    %find elements within a small region
    [enrDomain] = crackDetect(xCrk,ipas,tipElem,splitElem,vertexElem,cornerElem,enrDomain) ;

    %find type of element: tip, split, vertex
    [typeElem,elemCrk,tipElem,splitElem,vertexElem,cornerElem,tangentElem,xTip,xVertex,enrichNode,crackNode] = nnodeDetect(xCrk,enrDomain) ;
    % Deal with corner nodes by introducing phantom nodes (1 for each signed distance)
    
      %warning('Phantom nodes introduced to account for crack going through nodes')
      %[enrichNode, n_red] =  f_phantomNode(crackNode,elemCrk,splitElem,tipElem,vertexElem,enrichNode) ;
    %else
      %n_red = 0;
    %end 


    


        %plot enriched nodes
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

    [F] = ForceVector(F) ;

    else

      [K,F] = KmatXFEM3(enrichNode,elemCrk,typeElem,xTip,xVertex,...
        splitElem,tipElem,vertexElem,cornerElem,crackNode,enrDomain,pos,xCrk,K,F) ;

    end


    
    if exist('rift_wall_pressure') & strcmp(rift_wall_pressure,'y')
      [F] = f_apply_ocean_pressure(enrichNode,elemCrk,typeElem,xTip,xVertex,...
        splitElem,tipElem,vertexElem,cornerElem,crackNode,enrDomain,pos,xCrk,F) ;
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

    if any(fixedF)
      [crackLips,flagP] = f_cracklips( zeros(totalUnknown,1), xCr, enrDomain, typeElem, elemCrk, xTip,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);
      Fcrack = zeros(size(F));
      Fcrack = f_crackforce_fixed(Fcrack,fixedF,crackLips,xCr,elemCrk,xTip,pos,typeElem, enrichNode,splitElem,vertexElem,tipElem);
      F = F + Fcrack;
    end

    [L,U] = lu(K) ;
    y = L\F;
    u = U\y;

    fu = full(u);
    Stdux = fu(1:2:2*numnode) ;
    Stduy = fu(2:2:2*numnode) ;
    [crackLips,flagP] = f_cracklips( u, xCr, enrDomain, typeElem, elemCrk, xTip,xVertex,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);

    
    f = figure('visible','on');
    hold on
    dfac = 1 ;
    plotMesh(node+dfac*[Stdux, Stduy],element,elemType,'b-',plotNode,f)
    f_plotCrack(crackLips,1,'r-','g-','k--')
    print([results_path,'/crack_walls',num2str(ipas)],'-dpng','-r300')



    if contact & ~flagP
      % first we need to find out if there is any interpenetration
      penalty = 0
      disp([num2str(toc),'    No contact therefore penalty method was not applied'])
    end
    f = figure('visible','on');
    clf
    trisurf(element,node(:,1),node(:,2),Stduy)
    axis equal; view(2); shading interp; colorbar
    title('Y displacement before Newton solver')
    print('original_disp','-dpng')
    %keyboard




    if penalty
      tol = 1e-8;
      cont = 1
      Du = zeros(size(u));
      nu = 1
      while 1 
        disp(['Newton step ',num2str(cont)])
        disp(['---------------------------------------------'])
        Fint = K*u;
        Fext = F;
        [KT,Gint] = KTmatXFEM(Kpen,enrichNode,elemCrk,typeElem,xTip,xVertex,splitElem,tipElem,vertexElem,cornerElem,crackNode,pos,xCrk,K,u);
        Res  = Fext - Fint - Gint;
        nr = norm(Res,2);
        if cont == 1
          nr0  = nr;
        end
        rnr = nr/nr0;
        disp(['L2 norm of the residual, R =  ',num2str(nr)])
        disp(['Relative to R0 : ',num2str(rnr)])

        [L,U] = lu(KT) ;
        y = L\Res;
        Du = U\y;
        u = u + Du;
        cont = cont + 1;
        fu = full(u);
        Stdux = fu(1:2:2*numnode) ;
        Stduy = fu(2:2:2*numnode) ;

        [crackLips,flagP] = f_cracklips( u, xCr, enrDomain, typeElem, elemCrk, xTip,enrichNode,crackNode,pos,splitElem, vertexElem, tipElem);

        f = figure('visible','on');
        clf
        hold on
        dfac = 1 ;
        plotMesh(node+dfac*[Stdux, Stduy],element,elemType,'b-',plotNode,f)
        f_plotCrack(crackLips,1,'r-','g-','k--')
        print(['crack_iter',num2str(cont)],'-dpng','-r300')
        %keyboard


        if rnr < tol
           disp(['Converged at step : ',num2str(cont)])
           break
        elseif cont > 10
           warning(['After, ',num2str(cont),' iterations ||R||/||R0|| is still: ',num2str(rnr)])
           break
        end
      end
    end
%     
%     % plot displacement contour
     %figure
     %clf
     %trisurf(element,node(:,1),node(:,2),Stduy)
     %axis equal; view(2); shading interp; colorbar
     %title(['Y displacement after Newton solver (',num2str(cont) ,' iterations)'])
     %keyboard
     
     %save('test.mat','K','F','u')


    
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

    
     f = figure();
     hold on
     dfac = 50 ;
     plotMesh(node+dfac*[Stdux, Stduy],element,elemType,'b-',plotNode,f)
      %plotMesh(node+dfac*[uxAna uyAna],element,elemType,'r-',plotNode)
     figure_name = ['Disp_fact_',num2str(dfac),'_',num2str(ipas)];
     print([results_path,'/',figure_name],'-dpng','-r300')


   
    [Knum,Theta,xCrk] = KcalJint(xCrk,...
        typeElem,enrDomain,elemCrk,enrichNode,crackNode,xVertex,...
        vertexElem,pos,u,ipas,delta_inc,Knum,Theta,tipElem,splitElem,cornerElem) ;


end
save([results_path,'/crack_disp.mat'],'xCrk','Knum','Theta','u','element','node');
