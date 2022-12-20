function [sigma] = f_getstress(iel);
global  node element max_elem
global ISSM_xx ISSM_yy ISSM_xy

if ~isempty(ISSM_xx) 
  xx = ISSM_xx(iel); 
  yy = ISSM_yy(iel); 
  xy = ISSM_xy(iel); 
  sigma = [ 2*xx+yy xy; xy 2*yy+xx ];
else
  % get the centre of each element
  nodes = element(iel,:);
  xc    = mean(node(nodes,:),1);
  %extract stress from the stress field being applied 
  sigma = f_extractStress(xc)';
end

end

