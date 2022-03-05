function [type_elem,elem_crk,tip_elem,split_elem,vertex_elem,corner_elem...
     xTip,xVertex,enrich_node,crack_node] = nodeDetect(xCr,elems)

global node element epsilon

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

%f = figure()
%hold on
%plotMesh(node,element,'T3','b-','no')


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
        flag1 = 0;
        flag2 = 0;
        flag3 = 0;
        flag4 = 0;
        flag5 = 0;
        for kj = 1:size(xCr(kk).coor,1)-1       %loop over the elements of the fracture
            int_points = [];
            q1 = xCr(kk).coor(kj,:); 
            q2 = xCr(kk).coor(kj+1,:);
            for iedge=1:size(sctr,2)             %loop over the edges of elements
                nnode1=sctrl(iedge);
                nnode2=sctrl(iedge+1);
                p1 = [node(nnode1,:)];
                p2 = [node(nnode2,:)];
                intersect =segments_int_2d(p1,p2,q1,q2) ;
                intes = intes + intersect(1);
                if intersect(1)                
                  int_points = [ int_points ; intersect(2:3) ];
                end
              end
              
              if ~isempty(int_points) 
                flag1 = inhull(q1,vv,[],-1*epsilon);
                flag2 = inhull(q2,vv,[],-1*epsilon);
                [crk_int,flag3,crack_node] = f_update_crk_int(crk_int,int_points,flag1,flag2,q1,q2,crack_node,sctr);
                if kj == 1 & flag1
                  flag4 = 1;
                end
                if (kj == size(xCr(kk).coor,1)-1) & flag2
                  flag5 = 1;
                end
                xCr_element(e,:) = xCr(kk).coor(kj,:) * flag1 + xCr(kk).coor(kj+1,:) * flag2;  % link between crack coordinate and elements  
              end
              if iel == 3401
                keyboard
              end
        end
        

        %figure(f)
        %pvv = [vv;vv(1,:)] 
        %celem = plot(pvv(:,1),pvv(:,2),'r--')
                    %keyboard
        %delete(celem)
        
        % this is now true for all conditions
   %------- let's choose the categorie -------%     
       if (flag3 == 1)     % CORNER or SPLITCORNER or VERTEX CORNER 
         elem_crk(e,:) = crk_int;
         corner_elem = [corner_elem; e];
       end
       cond1 = ((intes == 2 ) & flag3 ~= 1 );
       cond2 = (intes == 3 );
       if ( cond1 | cond2 ) & (flag1 == 0) & (flag2 == 0)      % SPLIT  
            type_elem(elems(iel),kk) = 2; 
            split_elem = [split_elem; e];
            elem_crk(e,:) = crk_int;
        end
        if (((flag1 == 1) | flag2==1) & (intes>=2)) & not(flag4 | flag5)        % VERTEX      
            type_elem(e,kk) = 3; 
            vertex_elem = [vertex_elem; e] ;
            elem_crk(e,:) = crk_int;
            xVertex(e,:) = xCr_element(e,:); 
        end
        if  ( flag4 | flag5 )                                    % TIP
            type_elem(e,kk) = 1;  
            tip_elem = [tip_elem; e] ;
            xTip(e,:) = xCr_element(e,:);            
            elem_crk(e,:) = crk_int;  % coordinates needed for SIF computation           
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
