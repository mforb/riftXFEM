function [xCr] = f_move_vertex(xCr,move_kj);
persistent direction
global epsilon

if isempty(direction)
q1 = xCr.coor(move_kj,:)
q2 = xCr.coor(move_kj+1,:)
q3 = xCr.coor(move_kj+2,:)

seg1 = [q1,q2];
seg2 = [q2,q3];

v1 = seg1 - q1;
v2 = seg2 - q2;

n1 = norm(v1);
n2 = norm(v2);

v3 = n1*(v2/n2) + n2(v1/n1);
n3 = norm(v3);

direction = v3/n3;
end

xCr.coor(move_kj+1,:) = xCr.coor(move_kj+1,:) + direction*2*epsilon; %that simple!


