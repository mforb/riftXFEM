function pos = posi(xCr,numnode,enrich_node,crack_nodes)

global node element
global orig_nn

pos = zeros(numnode,size(xCr,2));
nsnode = 0 ;
ntnode = 0 ;

for k = 1:size(xCr,2)
    for i = 1 : numnode
        if (enrich_node(i,k) == 2)  %split element cuz of crack
            pos(i,k) = (numnode + nsnode*1 + ntnode*4) + 1 ;
            nsnode = nsnode + 1 ; % because otherwise we already have added degrees of freedom
            %if ~ismember(i,crack_nodes)
              %pos(i,k) = (numnode + nsnode*1 + ntnode*4) + 1 ;
              %nsnode = nsnode + 1 ; % because otherwise we already have added degrees of freedom
            %else
              %id = find(crack_nodes==i);
              %pos(i,k) = orig_nn + id
            %end
        elseif(enrich_node(i,k) == 3)   %split cuz of material
            pos(i,k) = (numnode + nsnode*1 + ntnode*4) + 1 ;
            nsnode = nsnode + 1 ;
        elseif (enrich_node(i,k) == 1)  %tip cuz of crack
            pos(i,k) = (numnode + nsnode*1 + ntnode*4) + 1 ;
            ntnode = ntnode + 1 ;
        end
    end
end
