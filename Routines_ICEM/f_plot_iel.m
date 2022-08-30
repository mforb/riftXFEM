function [ ] = f_plot_iel( iel, se )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Mon AUG 29 2022

global node element elemType

if ( nargin < 2 )
   se='b-';
end

x = node(element(iel,:),1);
y = node(element(iel,:),2);
x = [x;x(1)];
y = [y;y(1)];

plot(x,y,se)

