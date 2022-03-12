function [F] = f_getHeight(iel);
global  node element max_elem
global ISSM_xx ISSM_yy ISSM_xy
global rhow rhoi 

if ~isempty(ISSM_h) 
  H = [ISSM_xx(iel), ISSM_yy(iel), ISSM_xy(iel)]';
else
  % get the centre of each element
  nodes = element(iel,:);
  xc    = mean(node(nodes,:),1);
  %extract stress from the stress field being applied 
  H = f_extractHeight(xc)';
end

if ~exist('rhow') | ~exist('rhoi')
  rw = 1027;
  ri = 917;
else
  rw = rhow;
  ri = rhoi;
end

F = H*H*rhoi*(rhow-rhoi)/rhow

end

