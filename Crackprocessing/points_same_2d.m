function flag = points_same_2d( p1, p2 )
global epsilon
% points_same_2D tests if 2 points in 2D are the same
%
%  Parameters:
%
%    Input, real P1(2) 
%    Input, real P2(2)
%
%    Output, integer FLAG, records the results.
%    0, the points are not the same 
%    1, the points are the same 
tol = epsilon;
flag = 1;
u = abs(p1(1)-p2(1));
v = abs(p1(2)-p2(2));
if ( u > tol | v > tol )
  flag = 0;
end
end
