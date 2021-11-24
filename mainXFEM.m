function [Knum,Theta,xCrk] = mainXFEM(xCrk,npas,delta_inc)

%-- Declare global variables here
global elemType stressState typeCrack
global L D E nu C P sigmato
global numcrack xCr deltaInc numstep
global plotmesh plotNode
global node element numnode numelem bcNodes edgNodes typeProblem

Knum = [ ] ; Theta = [ ] ;
Enrdomain = [ ] ; tipElem = [ ] ; splitElem = [ ] ; vertexElem = [ ] ; cornerElem = [];
%loop over number of steps of crack growth
for ipas = 1:npas
    disp([num2str(toc),'    Crack growth number     ',num2str(ipas)]) ;

    disp([num2str(toc),'    Crack processing']) ;
    %find elements within a small region
    [Enrdomain] = crackDetect(xCrk,ipas,tipElem,splitElem,vertexElem,cornerElem,Enrdomain) ;

    %find type of element: tip, split, vertex
    [typeElem,elemCrk,tipElem,splitElem,vertexElem,cornerElem,xTip,xVertex,enrichNode,crackNode]=...
        nodeDetect(xCrk,Enrdomain) ;
    % Deal with corner nodes by introducing phantom nodes (1 for each signed distance)
    if ~isempty(crackNode)
      warning('Phantom nodes introduced to account for crack going through nodes')
      [enrichNode, n_red] =  f_phantomNode(crackNode,elemCrk,splitElem,tipElem,vertexElem,enrichNode) ;
    else
      n_red = 0;
    end 


    


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
    totalUnknown = numnode*ndof + split*1*ndof + tip*4*ndof - n_red*1*ndof ;

    K = sparse(totalUnknown,totalUnknown) ;
    F = sparse(totalUnknown,1) ;


    %stiffness matrix computation
    pos = posi(xCrk,numnode,enrichNode,crackNode) ;

    disp([num2str(toc),'    Stiffness Matrix Computation']) ;

    [K] = KmatXFEM(enrichNode,elemCrk,typeElem,xTip,xVertex,...
        splitElem,tipElem,vertexElem,cornerElem,crackNode,pos,xCrk,K) ;

    [F] = ForceVector(F) ;

    kk1 = K ;
    %----- Imposing Essential boundary conditions
    disp([num2str(toc),'    Imposing Essential boundary conditions'])
    
    if strcmp(typeProblem,'eCrkTen')
        
        dispNodes = unique([bcNodes{1}]) ;
        bcdof = [ ]; bcval = [ ];
        for i=1:length(dispNodes)
            bcdof = [bcdof 2*dispNodes(i)] ;
            bcval = [bcval 0] ;
        end
        bcdof = [bcdof 2*dispNodes(1)-1] ;
        bcval = [bcval 0];
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

    [L,U] = lu(K) ;
    y = L\F;
    u = U\y;

    u = full(u);
    Stdux = u(1:2:2*numnode) ;
    Stduy = u(2:2:2*numnode) ;
%     
%     % plot displacement contour
     figure
     clf
     trisurf(element,node(:,1),node(:,2),Stduy)
     axis equal; view(2); shading interp; colorbar
     title('Displacement from XFEM solution')
     
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
    
    
    
     figure
     hold on
     dfac = 50 ;
     plotMesh(node+dfac*[Stdux, Stduy],element,elemType,'b-',plotNode)
     %plotMesh(node+dfac*[uxAna uyAna],element,elemType,'r-',plotNode)
   
    [Knum,Theta,xCrk] = KcalJint(xCrk,...
        typeElem,Enrdomain,elemCrk,enrichNode,crackNode,xVertex,...
        vertexElem,pos,u,ipas,delta_inc,Knum,Theta,tipElem,splitElem,cornerElem) ;

end
