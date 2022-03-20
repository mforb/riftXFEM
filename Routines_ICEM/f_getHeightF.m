function [F] = f_getHeightF(iel);
global  node element max_elem
global ISSM_xx ISSM_yy ISSM_xy
global rhow rhoi 

g = 9.81
if exist('ISSM_h') 
  H = ISSM_h(iel);
else
  % get the centre of each element
  nodes = element(iel,:);
  xc    = mean(node(nodes,:),1);
  %extract stress from the stress field being applied 
  H = f_extractHeight(xc)';
end

if isempty(rhow) | isempty(rhoi)
  rw = 1027;
  ri = 917;
else
  rw = rhow;
  ri = rhoi;
end

F = -1*H*g*ri*(rw-ri)/rw;

end

