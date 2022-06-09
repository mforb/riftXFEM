function [type_elem,elem_crk,tip_elem,split_elem,vertex_elem,corner_elem,tangent_elem...
     xTip,xVertex,enrich_node,crack_node,xCr] = nnodeDetect(xCr,elems)

global node element epsilon
global plothelp elemType
global output_file

type_elem = zeros(size(element,1),size(xCr,2)) ;
elem_crk = zeros(size(element,1),4) ;
xTip = zeros(size(element,1),2) ;
xVertex = zeros(size(element,1),2) ;
enrich_node = zeros(size(node,1),size(xCr,2)) ;
crack_node  = [];
tip_elem = [] ;
split_elem = [] ;
vertex_elem = [] ;
corner_elem = [] ; 
tangent_elem = [] ; 
crack_node = [] ; % this can be improved for a more efficient step by step code

%hold on
%plotMesh(node,element,'T3','b-','no')
if plothelp
  f = figure()
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
  mv_kj_old = 0;
  clear f_move_vertex;
  while 1
    mv_kj = 0;
    crack_node  = [];
    corner_elem = [] ; 
    tangent_elem = [] ; 
    subelems = [];
    for i = 1:length(elems)
      [subelems, tangent_elem, corner_elem, type_elem, elem_crk, crack_node, kj_track ] = f_find_subelems(elems(i),xCr(kk),kk,subelems, tangent_elem, corner_elem, type_elem, elem_crk, crack_node);
      if size(kj_track,2)>1
        warn_str = ['Element ',num2str(elems(i)),': interacts with mesh twice at intersection between crack segmenets ',num2str(kj_track(1)),' and ',num2str(kj_track(2))];
        warning(join(warn_str));
        fprintf(output_file,[warn_str,'\n'])
        mv_kj = kj_track(1);
        break;
      end
    end

    if mv_kj == 0
      % success!
      break;
    else
      if mv_kj ~= mv_kj_old
        clear f_move_vertex % we want to find a new direction to move vertex kj
      end
      [xCr_new] = f_move_vertex(xCr(kk),mv_kj);
      mv_kj_old = mv_kj;
      xCr(kk) = xCr_new;
      % this can be improved if necessary
      continue;
    end
  end

  found_start = 0;
  found_end   = 0;

  for iel=1:length(subelems)                     %loop on elems (=elements selected for enrichment)
    e = subelems(iel) ;
    sctr=element(e,:);
    sctrl = [sctr sctr(1,1)];      
    vv = node(sctr,:);
    crk_int = [];
    flag1 = 0; % is the first point of a crack segment is in this element
    flag2 = 0; % is last point of the crack segment is in the element
    flag3 = 0; % there is a first or last crack point in the element 
    flag4 = 0; % is there one or more nodes that are coincident with the crack (tangent to the crack)
    flag5 = 1; % split element 
    for kj = 1:size(xCr(kk).coor,1)-1       %loop over the elements of the fracture
      q1 = xCr(kk).coor(kj,:); 
      q2 = xCr(kk).coor(kj+1,:);
      if ~flag1
        flag1 = inhull(q1,vv,[]) ; % we need to detect these in order to run the simulation
        if flag1
          if (kj == 1 ) & flag1 
            flag3 = flag1 ;
          end
          flag5 = 0;
          xCr_element = q1;  % link between crack coordinate and elements  
        end
      end
      if ~flag2
        flag2 = inhull(q2,vv,[]) ; % otherwise we might not have end points
        if flag2
          if (kj == size(xCr(kk).coor,1)-1) & flag2
            flag4 = flag2;
          end
          flag5 = 0;
          xCr_element = q2;
        end
      end

            % find the intersection points
      int_points = [];

      for iedge=1:size(sctr,2)             %loop over the edges of elements
          nnode1=sctrl(iedge);
          nnode2=sctrl(iedge+1);
          p1 = [node(nnode1,:)];
          p2 = [node(nnode2,:)];
          intersect = segments_int_2d(p1,p2,q1,q2) ;
          %edge_track=[edge_track,iedge]
          if intersect(1)                
            int_points = [ int_points ; intersect(2:3) ];
          end
      end
              
      if ~isempty(int_points)
        [crk_int] = f_up_crk_int(crk_int,int_points,flag1,flag2,q1,q2,e);
      end
    end % crack segment
    if flag5 & isempty(crk_int)
      %disp(['Should have picked up an intersection, check element ',num2str(e)])
      flag5 = 0; %Can be an element around the tip element
    end

    %figure(f)
    %pvv = [vv;vv(1,:)] 
    %celem = plot(pvv(:,1),pvv(:,2),'r--')
    %if flag1 | flag2
      %keyboard
    %end
    %delete(celem)
    
    % this is now true for all conditions
 %------- let's choose the categorie -------%   
     % check for tangent elements
    if flag5     % SPLIT  
      type_elem(e,kk) = 2; 
      split_elem = [split_elem; e];
      elem_crk(e,:) = crk_int;
    elseif flag1 & flag2         % VERTEX      
      type_elem(e,kk) = 3; 
      vertex_elem = [vertex_elem; e] ;
      xVertex(e,:) = xCr_element; 
      elem_crk(e,:) = crk_int;
    elseif   flag3                                     % TIP
      found_start = 1;
      type_elem(e,kk) = 1;  
      tip_elem = [tip_elem; e] ;
      xTip(e,:) = xCr_element;            
      elem_crk(e,:) = crk_int;
    elseif flag4
      found_end = 1;
      type_elem(e,kk) = 1;  
      tip_elem = [tip_elem; e] ;
      xTip(e,:) = xCr_element;            
      elem_crk(e,:) = crk_int;
    end

    if ismember(e,corner_elem)
      for ni = 1:length(sctr) 
        no = sctr(ni);
        if ismember(no,crack_node)
          try
            d1 = norm(crk_int(1:2)-node(no,:),2);
            d2 = norm(crk_int(3:4)-node(no,:),2);
          catch
            keyboard
          end
          if d1<d2
            elem_crk(e,1:2) = node(no,:);
          else
            elem_crk(e,3:4) = node(no,:);
          end
        end
      end
    end

  end % iel    
end % kk

if ~found_start
  warning('No crack starting point found')
end
if ~found_end
  warning('No crack ending point found')
end

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
                if enrich_node(sctr(in),kk) == 1 
                  warning(['tangent element next to tip element (might not work) at:  ',num2str(iel)])
                elseif enrich_node(sctr(in),kk) == 0 
                            enrich_node(sctr(in),kk)   = 2;
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

