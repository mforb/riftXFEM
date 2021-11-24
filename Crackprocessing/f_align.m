function [ pp ] = f_align( p1, p2, xdir )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Fri Nov 12 19:13:22 NZDT 2021
pp1 = p1 - p2;
dpp1 = pp1(1)*xdir(1)+pp1(2)*xdir(2);
if dpp1 > 0
  pp = [ p2,p1];
else
  pp = [ p1, p2 ];
end
end
