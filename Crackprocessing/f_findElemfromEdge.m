function [ elem ] = f_findElemfromEdge( edge )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Tue Apr 12 09:26:59 NZST 2022
global node element

% the edge is composed of two nodes, there is only one element that will contain both those nodes
[e1,r] = find(element == edge(1));
[e2,r] = find(element == edge(2));
elem = intersect(e1,e2);
end
