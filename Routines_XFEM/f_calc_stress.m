function [sigma] = f_calc_stress(pt,e,u,C,type_elem,enrich_node,elem_crk,xVertex,xyTip,crack_nodes,pos,k,varargin)

global node element elemType E nu Cm1

if nargin == 13
  QT = varargin{1};
else
  QT = eye(2);
end
B = [];
B = [B xfemBmat(pt,e,type_elem,enrich_node(:,k),elem_crk,xVertex,xyTip,crack_nodes,k)];

leB = size(B,2);

U = [];
U     = [U; element_disp(e,pos(:,k),enrich_node(:,k),u,k)];

epsilon = B*U ;
sigma   = C*epsilon;

voit2ind    = [1 3;3 2];
sigma   = QT*sigma(voit2ind)*QT';
