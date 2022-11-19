function [xCr,theta_inc,kstr] = f_advance_tip(xCr,tip,delta_inc,theta_inc,alpha,kstr)
persistent ang1
persistent ang2
if tip == 1
  if isempty(ang1)
    ang1 = 0.01;
  else
  ang1 = -ang1;
  end
 ang = ang1
 cd = 1;
 c_fld = 'coornew1';
elseif tip == 2
  if isempty(ang2)
    ang2 = 0.01;
  else
    ang2 = -ang2;
  end
 ang = ang2
 cd = size(xCr.coor,1) ;
 c_fld = 'coornew2';
end


inc_x = xCr.coor(cd,1) + delta_inc * (cos(theta_inc)*cos(alpha) - sin(theta_inc)*sin(alpha));
inc_y = xCr.coor(cd,2) + delta_inc * (cos(theta_inc)*sin(alpha) + sin(theta_inc)*cos(alpha)); % y is flipped in the coordinate system used in SIF so that the "positive" side is the same
tt = f_test_tip([inc_x,inc_y],delta_inc)
while tt
  theta_inc = theta_inc + ang;
  inc_x = xCr.coor(cd,1) + delta_inc * (cos(theta_inc)*cos(alpha) - sin(theta_inc)*sin(alpha));
  inc_y = xCr.coor(cd,2) + delta_inc * (cos(theta_inc)*sin(alpha) + sin(theta_inc)*cos(alpha)); % y is flipped in the coordinate system used in SIF so that the "positive" side is the same
  tt = f_test_tip([inc_x,inc_y],delta_inc)
  kstr = [kstr, 'Tip ',num2str(tip),': Theta was modified by ',num2str(ang),' to avoid tip coincident with a node\n']; 
end

xCr.(c_fld)= [inc_x inc_y]; %
