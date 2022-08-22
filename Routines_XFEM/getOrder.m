function [intOrder]=getOrder(iel,splitElem,tipElem,vertexElem,tipEnr,tipNodes,...
    TipOrder,SplitOrder,VertexOrder,NormalOrder)

%declare global variables here
global node element

sctr = element(iel,:) ;
nnode = node(sctr,:) ;

if( ismember(iel,tipElem) )      %tip element
    intOrder = TipOrder ;
elseif( ismember(iel,splitElem) & any(ismember(sctr,tipNodes)) )
    %split element with tip nodes
    intOrder = TipOrder ;
elseif( ismember(iel,vertexElem) )       %vertex element
    intOrder = VertexOrder ;
    subTriDiv = 0 ;
elseif( size(tipEnr,1) > 0 )
    %standard Q4 element with tip sharing nodes
    intOrder = TipOrder ;
elseif( ismember(iel,splitElem) )
    %split element..... only heaviside enrichment
    for c = 1:length(splitElem)
        M = find(iel == splitElem(c) ) ;
        if( M == 1)
            intOrder = SplitOrder ;
        end
    end
else        %Standard Q4 element
   intOrder = NormalOrder ;
end