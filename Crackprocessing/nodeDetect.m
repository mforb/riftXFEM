function [type_elem,elem_crk,tip_elem,split_elem,vertex_elem...
    ,xTip,xVertex,enrich_node] = nodeDetect(xCr,elems)

global node element

type_elem = zeros(size(element,1),size(xCr,2)) ;
elem_crk = zeros(size(element,1),4) ;
xCr_element = zeros(size(element,1),2) ;
xTip = zeros(size(element,1),2) ;
xVertex = zeros(size(element,1),2) ;
enrich_node = zeros(size(node,1),size(xCr,2)) ;
tip_elem = [] ;
split_elem = [] ;
vertex_elem = [] ;

%select the special elements(tip, vertex, split)    %loop on
%elems(=elements selected for enrichment)
% select the special elements(tip, vertex, split)
for kk = 1:size(xCr,2)  
    for iel=1:size(elems,1)                     %loop on elems (=elements selected for enrichment)
        e = elems(iel) ;
        sctr=element(e,:);
        vv = node(sctr,:);
        crk_int = [];
        intes = 0;
        flag1 = 0;
        flag2 = 0;
        for kj = 1:size(xCr(kk).coor,1)-1       %loop over the elements of the fracture
            q1 = xCr(kk).coor(kj,:); 
            q2 = xCr(kk).coor(kj+1,:);
            sctrl = [sctr sctr(1,1)];      
            for iedge=1:size(sctr,2)             %loop over the edges of elements
                nnode1=sctrl(iedge);
                nnode2=sctrl(iedge+1);
                p1 = [node(nnode1,:)];
                p2 = [node(nnode2,:)];
                intersect=segments_int_2d(p1,p2,q1,q2) ;
                intes = intes + intersect(1);
                if intersect(1) > 0
                    crk_int = [crk_int intersect(2) intersect(3)];
                    flag1 = inhull(xCr(kk).coor(kj,:),vv,[],-1e-8);
                    flag2 = inhull(xCr(kk).coor(kj+1,:),vv,[],-1e-8);
                    xCr_element(e,:) = xCr(kk).coor(kj,:) * flag1 + xCr(kk).coor(kj+1,:) * flag2;  % link between crack coordinate and elements  
                end
            end %for iedge
        end
        
   %------- let's choose the categorie -------%     
        if ((intes == 2) & (flag1 == 0) & (flag2 == 0))     % SPLIT
            type_elem(elems(iel),kk) = 2; 
            split_elem = [split_elem; e];
            elem_crk(e,:) = crk_int;
        end
        if (((flag1 == 1) | flag2==1) & (intes==2))         % VERTEX      
            type_elem(e,kk) = 3; 
            vertex_elem = [vertex_elem; e] ;
            elem_crk(e,:) = crk_int;
            xVertex(e,:) = xCr_element(e,:); 
        end
        if (intes == 1)                                     % TIP
            type_elem(e,kk) = 1;  
            tip_elem = [tip_elem; e] ;
            xTip(e,:) = xCr_element(e,:);
            elem_crk(e,:) = [crk_int xTip(e,1) xTip(e,2)];  % coordinates needed for SIF computation           
        end
    end % iel    
end % kk

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
        end  %if
    end  %loop on the elements
end  %loop on cracks
