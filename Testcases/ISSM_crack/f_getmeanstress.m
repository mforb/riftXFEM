function [m_stress] = f_getmeanstress(element);
sigma = [];
for iel = 1:size(element,1)
  sigma = [ sigma; getstress(iel)']; 
end

m_stress = mean(sigma,1);
