function [ fd_xy, fu_xy ] = f_calc_crack_force( p , cd , cu, F_app )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Mon Dec  6 19:35:27 NZDT 2021
  [W,Q] = quadrature(1,'GAUSS',1) ;
  p1 = p(1:2);
  p2 = p(3:4);
  ddown1 = cd(1:2);
  ddown2 = cd(3:4);
  dup1   = cu(1:2);
  dup2   = cu(3:4);
  c_d1 = p1 + ddown1;  
  c_p1 = p1 + dup1;
  c_d2 = p2 + ddown2;
  c_p2 = p2 + dup2;
  c_dtan = c_d2 - c_d1;
  c_ptan = c_p2 - c_p1;
  c_dnorm = sqrt(c_dtan(1)^2 +c_dtan(2)^2);  
  c_pnorm = sqrt(c_ptan(1)^2 +c_ptan(2)^2);  
  c_dm   = c_d1 + 0.5*c_dtan;
  c_pm   = c_p1 + 0.5*c_ptan;
  c_dtan = c_dtan/c_dnorm;  
  c_ptan = c_ptan/c_pnorm;  
  % first we calculate the normal and tangental force on each corner of the crack lips
  fu = 0;
  fd = 0;
  fut = 0;
  fdt = 0;
  for q = 1:size(W,1)
    pt = Q(q,:) ; wt = W(q) ;
    N = lagrange_basis('L2',pt) ;
    J_down = c_dnorm/2 ;
    J_up = c_pnorm/2 ;
    fu  = fu + N*F_app(2)*J_up*wt ;
    fd  = fd + N*F_app(2)*J_down*wt ;
    fut  = fut + N*F_app(1)*J_up*wt ;
    fdt  = fdt + N*F_app(1)*J_down*wt ;
  end
  % Now we figure out what the forces are in x and y
  alpha1 = atan2(c_dtan(2),c_dtan(1));
  alpha2 = atan2(c_ptan(2),c_ptan(1));
  QT1  =[cos(alpha1) sin(alpha1); -sin(alpha1) cos(alpha1)];          
  QT2  =[cos(alpha2) sin(alpha2); -sin(alpha2) cos(alpha2)];          
  fu1 = [fut(1) fu(1)];
  fu2 = [fut(2) fu(2)];
  fd1 = [-fdt(1) -fd(1)]; %not planning on using tangential forces but care is need in choosing direction
  fd2 = [-fdt(2) -fd(2)];
  fu_xy = [ fu1*QT2 , fu2*QT2 ];
  fd_xy = [ fd1*QT1, fd2*QT1 ];
end
