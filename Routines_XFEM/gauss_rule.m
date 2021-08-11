function [W,Q] = gauss_rule(e,enrich_node,elem_crk,...
    xTip,xVertex,tip_elem,split_elem,vertex_elem,xCrk)

%get integration points and weights for elements

%declare global variables here
global node element

sctr = element(e,:) ;
nnode = node(sctr,:) ;

tip_enr = find(enrich_node(sctr,:) == 1) ;
tipnodes = [ ] ;
for k=1:size(xCrk,2)
    tipnodes = [tipnodes; find(enrich_node(:,k) == 1)];
end

[phi] = dista(e,elem_crk) ;

%numerical integration order for numerical integration
NormalOrder = 2 ;  %max = 8
TipOrder = 6 ;     %max = 20
SplitOrder = 3;
VertexOrder = 3 ;

[IntOrder]=getOrder(e,split_elem,tip_elem,vertex_elem,tip_enr,tipnodes,...
    TipOrder,SplitOrder,VertexOrder,NormalOrder) ;

if( IntOrder == TipOrder && TipOrder <= 7 )
    intType = 'DUNAVANT' ;        %options: TRIANGULAR or DUNAVANT
    subTriDiv = 0 ;
elseif( IntOrder == TipOrder && TipOrder > 7)
    intType = 'DUNAVANT' ;
    subTriDiv = 0 ;
end

if( ismember(e,tip_elem) )      %tip element
    [W,Q] = disTipQ4(IntOrder,phi,nnode,xTip(e,:),subTriDiv,intType) ;
elseif( ismember(e,split_elem) && any(ismember(sctr,tipnodes)) )
    %split element with tip nodes
    [W,Q] = disSplitQ4(IntOrder,phi,subTriDiv,intType) ;
elseif( ismember(e,vertex_elem) )       %vertex element
    subTriDiv = 0 ;
    intType = 'TRIANGULAR' ;        %options: TRIANGULAR or DUNAVANT
    [W,Q] = disTipQ4(IntOrder,phi,nnode,xVertex(e,:),subTriDiv,intType) ;
elseif( size(tip_enr,1) > 0 )
    %standard Q4 element with tip sharing nodes
    [W,Q] = disBlendQ4(IntOrder,subTriDiv,intType) ;
elseif( ismember(e,split_elem) )
    %split element..... only heaviside enrichment
    for c = 1:length(split_elem)
        M = find(e == split_elem(c) ) ;
        if( M == 1)
            subTriDiv = 0 ;
            intType = 'TRIANGULAR' ;        %options: TRIANGULAR or DUNAVANT
            [W,Q] = disSplitQ4(IntOrder,phi,subTriDiv,intType) ;
        end
    end
else        %Standard Q4 element
    intType = 'GAUSS' ;
    [W,Q] = quadrature(IntOrder,intType,2) ;
end