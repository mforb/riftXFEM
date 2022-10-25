function [xCr] = f_move_vertex(xCr,move_kj);
persistent direction
global epsilon
global output_file

if isempty(direction)
  q1 = xCr.coor(move_kj,:)
  q2 = xCr.coor(move_kj+1,:)
  q3 = xCr.coor(move_kj+2,:)

  v1 = q1 - q2; 
  v2 = q3 - q2; 

  n1 = norm(v1);
  n2 = norm(v2);

  v3 = v2/n2 + v1/n1;
  n3 = norm(v3);

  direction = v3/n3;
end
kj_str = ['Crack vertex ',num2str(move_kj),': moved ',num2str(2*epsilon),' units in direction [',num2str(direction(1)),',',num2str(direction(2)),']'];
fprintf(output_file,[kj_str,'\n'])

xCr.coor(move_kj+1,:) = xCr.coor(move_kj+1,:) + direction*2*epsilon; %that simple!


