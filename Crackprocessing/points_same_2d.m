function flag = points_same_2d( p1, p2, tol )
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
if nargin == 2
tol = epsilon;
end

flag = 1;
u = p1(1)-p2(1);
v = p1(2)-p2(2);
r = sqrt(u*u+v*v)
if ( r > tol )
  flag = 0;
end
end
