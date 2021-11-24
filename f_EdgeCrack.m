function [ K, varargout ] = f_EdgeCrack( xCr,elemType,L,D,ndiv,)
% This is the function version of EdgeCrack, used mostly for testing the different run options
% Due to the (over)use of global variables. The variables are modified in this function (but need to be created beforehand)
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Thu Aug 12 17:33:16 NZST 2021
nout = max(nargout,1) - 1;
if strcmp(elemType,'Q4')
  rd = 0.0 ;
  [node,element,bcNodes,edgNodes] = createmesh(ndiv,rd) ;
  numnode = size(node,1) ;
  numelem = size(element,1) ;
elseif strcmp(elemType,'T3')
  load('Data/Mesh/T3_egde_test.mat')
end



[K,ThetaInc,xCr] = mainXFEM(xCr,numstep,deltaInc) 

if nargout == 1
  varargout{1} = ThetaInc;
elseif nargout == 2
  varargout{1} = ThetaInc;
  varargout{2} = xCr;
end

end
