function [ elem_force ] = f_readjust_elemforce( split_elem, elem_force)
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Fri Mar 18 18:39:27 NZDT 2022

%declare global variables here
global node element numnode numelem elemType




%loop over elements
for ii=1:size(split_elem,1)
  iel = split_elem(ii);
  fh = f_getHeightF(iel);
  elem_force(iel,[1,3]) = elem_force(iel,[1,3]) + fh ;
end
