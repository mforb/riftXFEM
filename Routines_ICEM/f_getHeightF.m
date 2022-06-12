function [F] = f_getHeightF(iel);
global  node element max_elem
global ISSM_H
global rhow rhoi 

g = 9.81;
if ~isempty(ISSM_H) 
  H = ISSM_H(iel);
else
  % get the centre of each element
  nodes = element(iel,:);
  xc    = mean(node(nodes,:),1);
  %extract height 
  try
    H = f_extractHeight(xc)';
  catch
    keyboard
  end
end

if isempty(rhow) | isempty(rhoi)
  rw = 1027;
  ri = 917;
else
  rw = rhow;
  ri = rhoi;
end

F = -0.5*H*H*g*ri*(1-ri/rw);

end

