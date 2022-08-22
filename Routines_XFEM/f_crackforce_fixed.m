function [ Fc ] = f_crackforce_fixed( Fc,F_app, crack_lips, xCr, xCrl, xTip, pos,type_elem, enr_node, split_elem, vertex_elem, tip_elem )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Mon Nov 29 11:58:33 NZDT 2021
global node element elemType

elems = union(split_elem,vertex_elem);
elems = union(elems,tip_elem);

for kk = 1:size(xCr,2) %what's the crack?
  for ii=1:size(elems,1)
    iel = elems(ii) ;
    sctr=element(iel,:);
    nn = length(sctr);
    p1 = crack_lips(ii,1:2,4,kk);
    p2 = crack_lips(ii,3:4,4,kk);
    p3 = crack_lips(ii,5:6,4,kk);
    n1 = crack_lips(ii,1:2,1,kk);
    n2 = crack_lips(ii,3:4,1,kk);
    n3 = crack_lips(ii,5:6,1,kk);
    dup1 = crack_lips(ii,1:2,2,kk);
    dup2 = crack_lips(ii,3:4,2,kk);
    dup3 = crack_lips(ii,5:6,2,kk);
    ddown1 = crack_lips(ii,1:2,3,kk);
    ddown2 = crack_lips(ii,3:4,3,kk);
    ddown3 = crack_lips(ii,5:6,3,kk);

    [fd_xy,fu_xy] = f_calc_crack_force([p1,p2],[ddown1,ddown2], [dup1,dup2],F_app);
    Fc = f_apply_crack_force(Fc,fd_xy,fu_xy,iel,xCr,xCrl,xTip,pos,type_elem,enr_node,n1,n2);


    if ismember(iel,vertex_elem)
      [fd_xy,fu_xy] = f_calc_crack_force([p2,p3],[ddown2,ddown3], [dup2, dup3], F_app);
      Fc = f_apply_crack_force(Fc,fd_xy,fu_xy,n2,n3,iel,xCr,xCrl,xTip,pos,type_elem,enr_node, n1, n2 );

    end
  end
end
