function [type_elem,elem_crk,tip_elem,split_elem,vertex_elem,corner_elem,tangent_elem...
     xTip,xVertex,enrich_node,crack_node] = nnodeDetect(xCr,elems)

global node element epsilon
global plothelp elemType

type_elem = zeros(size(element,1),size(xCr,2)) ;
elem_crk = zeros(size(element,1),4) ;
xCr_element = zeros(size(element,1),2) ;
xTip = zeros(size(element,1),2) ;
xVertex = zeros(size(element,1),2) ;
enrich_node = zeros(size(node,1),size(xCr,2)) ;
crack_node  = [];
tip_elem = [] ;
split_elem = [] ;
vertex_elem = [] ;
corner_elem = [] ; 
tangent_elem = [] ; 


%f = figure()
%hold on
%plotMesh(node,element,'T3','b-','no')
if plothelp
  figure(1)
  clf
  hold on
  plotMesh_numbered(node,element,elemType,'b-','no')
  for k=1:size(xCr,2)
      for kj = 1:size(xCr(k).coor,1)-1
          cr = plot(xCr(k).coor(kj:kj+1,1),xCr(k).coor(kj:kj+1,2),'k-') ;
          set(cr,'LineWidth',1);
      end
  end
end


%select the special elements(tip, vertex, split)    %loop on
%elems(=elements selected for enrichment)
% select the special elements(tip, vertex, split)
for kk = 1:size(xCr,2)  
    for iel=1:size(elems,1)                     %loop on elems (=elements selected for enrichment)
        e = elems(iel) ;
        sctr=element(e,:);
        sctrl = [sctr sctr(1,1)];      
        vv = node(sctr,:);
        crk_int = [];
        intes = 0;
        flag1 = 0; % is the first point of a crack segment is in this element
        flag2 = 0; % is last point of the crack segment is in the element
        flag3 = 0; % there is a first or last crack point in the element 
        flag4 = 0; % is there one or more nodes that are coincident with the crack (tangent to the crack)
        flag5 = 0; % does the crack intersect the element
        flag6 = 0; % tangent element
        for kj = 1:size(xCr(kk).coor,1)-1       %loop over the elements of the fracture
            q1 = xCr(kk).coor(kj,:); 
            q2 = xCr(kk).coor(kj+1,:);
            if epsilon > 1e-6
              in_epsi = 1e-6;
            else
              in_epsi = epsilon;
            end
            flag1 = inhull(q1,vv,[],-1*in_epsi) ; % we need to detect these in order to run the simulation
            flag2 = inhull(q2,vv,[],-1*in_epsi) ; % otherwise we might not have end points
            if flag1 | flag2
               xCr_element(e,:) = xCr(kk).coor(kj,:) * flag1 + xCr(kk).coor(kj+1,:) * flag2;  % link between crack coordinate and elements  
            end
            if (kj == 1 ) & flag1 
              flag3 = flag1 ;
            end
            if (kj == size(xCr(kk).coor,1)-1) & flag2
              flag3 = flag2;
            end

            %keyboard

            [ phi  ] = dista(e,[q1,q2]);
            [ psi1 ] = f_dista2(e,[q1,q2],q1);
            [ psi2 ] = f_dista2(e,[q2,q1],q2);
            lck_n = [];

            for nn = 1:length(vv)
                if abs(phi(nn)) < epsilon
                  if (abs(psi1(nn)) < epsilon) | (abs(psi2(nn)) < epsilon) | (psi1(nn)*psi2(nn) > 0)
                    crack_node = [ sctr(nn), crack_node];
                    flag4 = 1;
                    lck_n = [lck_n,nn];
                  end
                end
            end
            
            %keyboard

            % tangent elements
            if flag4
               nn = 1:length(vv)
               it_n = setdiff(nn,lck_n)
               if ~any(diff(sign(phi(it_n))))
                 flag6 = 1;
               end
            end

            if ~flag6
              % find the intersection points
              int_points = [];

              for iedge=1:size(sctr,2)             %loop over the edges of elements
                  nnode1=sctrl(iedge);
                  nnode2=sctrl(iedge+1);
                  p1 = [node(nnode1,:)];
                  p2 = [node(nnode2,:)];
                  intersect = segments_int_2d(p1,p2,q1,q2) ;
                  intes = intes + intersect(1);
                  %edge_track=[edge_track,iedge]
                  if intersect(1)                
                    int_points = [ int_points ; intersect(2:3) ];
                  end
              end
                
              if ~isempty(int_points)
                try 
                [crk_int] = f_up_crk_int(crk_int,int_points,flag1,flag2,q1,q2,e);
                catch
                  keyboard
                end
                flag5 = 1;
              end
            else
              crk_int = [q1,q2];
            end
        end
        

        

        %figure(f)
        %pvv = [vv;vv(1,:)] 
        %celem = plot(pvv(:,1),pvv(:,2),'r--')
                    %keyboard
        %delete(celem)
        
        % this is now true for all conditions
   %------- let's choose the categorie -------%     
       % check for tangent elements
       if flag4
          corner_elem = [corner_elem,e];
       end


       if flag6 
          type_elem(elems(iel),kk) = 4; 
          tangent_elem = [tangent_elem,e];
          elem_crk(e,:) = crk_int;
       elseif flag5 & ~(flag1 | flag2)    % SPLIT  
          type_elem(elems(iel),kk) = 2; 
          split_elem = [split_elem; e];
          elem_crk(e,:) = crk_int;
       elseif (flag1 | flag2 ) & ~flag3         % VERTEX      
          type_elem(e,kk) = 3; 
          vertex_elem = [vertex_elem; e] ;
          xVertex(e,:) = xCr_element(e,:); 
          elem_crk(e,:) = crk_int;
       elseif  ( flag3 )                                    % TIP
          type_elem(e,kk) = 1;  
          tip_elem = [tip_elem; e] ;
          xTip(e,:) = xCr_element(e,:);            
          elem_crk(e,:) = crk_int;
        end
    end % iel    
end % kk

crack_node = unique(crack_node);

% select the enriched nodes
for kk = 1:size(xCr,2)    
   for iel=1:size(elems,1)                     %loop on elems (=elements selected for enrichment)
        sctr = element(elems(iel),:);
        if type_elem(elems(iel),kk) == 1        % tip
            enrich_node(sctr,kk) = 1;
        elseif  type_elem(elems(iel),kk) == 2   % split
            for in=1:length(sctr)               % loop on the nodes of the element
                if enrich_node (sctr(in),kk) == 0  % already enriched
                    	enrich_node(sctr(in),kk)   = 2;
                end
            end
        elseif  type_elem(elems(iel),kk) == 4   % Tangent 
            for in=1:length(sctr)               % loop on the nodes of the element
                if enrich_node == 1 | enrich_node == 3
                  warning(['',num2str(iel)])

                end
                if enrich_node (sctr(in),kk) == 0  % already enriched
                    [Aw, Awp] = support_area(sctr(in),elems(iel),type_elem,elem_crk,xVertex,kk);
                    if (abs(Awp / Aw) > 1e-4) & (abs((Aw-Awp) / Aw) > 1e-4)   
                    	enrich_node(sctr(in),kk)   = 2;
                    end
                end
            end
        elseif type_elem(elems(iel),kk) == 3    %vertex 
            for in=1:length(sctr)
                if enrich_node (sctr(in),kk) == 0  % already enriched
                    [Aw, Awp] = support_area(sctr(in),elems(iel),type_elem,elem_crk,xVertex,kk); % let's test the tolerance linked to the support area! 
                    if ((abs(Awp / Aw) > 1e-4) & (abs((Aw-Awp) / Aw) > 1e-4)) 
                        enrich_node(sctr(in),kk)   = 2;
                    end
                end
            end  % loop on the nodes
        end
    end  %loop on the elements
end  %loop on cracks

