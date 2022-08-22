function [sigma] = f_getstress(iel);
global  node element max_elem
global ISSM_xx ISSM_yy ISSM_xy

if ~isempty(ISSM_xx) 
  sigma = [ISSM_xx(iel), ISSM_yy(iel), ISSM_xy(iel)]';
else
  % get the centre of each element
  nodes = element(iel,:);
  xc    = mean(node(nodes,:),1);
  %extract stress from the stress field being applied 
  sigma = f_extractStress(xc)';
end

end

