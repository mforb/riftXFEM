function p = f_align_lp_gc(p,elemck,sctr)
global node elemType

crkdir = elemck(3:4) - elemck(1:2);
[N1,dNdxi] = lagrange_basis(elemType,p(1,:)) ;
[N2,dNdxi] = lagrange_basis(elemType,p(end,:)) ;
pdir = N2' * node(sctr,:) - N1'*node(sctr,:);
if pdir*crkdir'< 0
  p([1,end],:) = p([end,1],:);
end
