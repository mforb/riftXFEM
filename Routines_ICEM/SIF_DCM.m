function [Knum,theta_inc] = SIF_DCM(C,flag_end,tip,elem_crk,xCr,type_elem,enrich_node,crack_nodes,xVertex,pos,u,F,kk,alpha,...
    tip_elem,split_elem,vertex_elem,corner_elem,tan_elem,elem_force)
    % main difference here is that F is also passed into the function

global node element elemType E nu Cm1
global iMethod iParam
global incR xc yc phiN
global lambda1 lambda2 nu1 nu2
global elemType typeMesh typeProblem typeCrack stressState
global wall_int output_file wall_force crack_load fixwf
global quick_freeze melange

plotNode = 'NO' ;


if strcmp(elemType,'Q4') 
  intType = 'GAUSS' ;
  corner = [1 2 3 4 1] ;
  nnode = [-1 -1;1 -1;1 1;-1 1] ;
else
  intType = 'TRIANGULAR';
  corner = [1 2 3 1] ;
  nnode = [0 0;1 0;0 1] ;
end
% for debug
%TR = triangulation(element,node);
%figure(12)
%triplot(TR)
%hold on
d = 1
% Compute the Stress Intensity Factors
% Using the Interaction integral method

% Steps :
% 1- detection of the elements on which we integrate
% 2- loop over these elements
% 3- loop over Gauss points
% 4- computation of stress, strain... in local coordinates !!!   ATTENTION
% 5- computation of the auxilliary fields: AuxStress and AuxEps and AuxGradDisp
% 6- computation of I1 and I2

% Determine J domain and weight function
if flag_end == 1
  xyTip = [elem_crk(tip,1) elem_crk(tip,2)] ;
  fl = [1 0; 0 -1];
elseif flag_end == 2
  xyTip = [elem_crk(tip,3) elem_crk(tip,4)] ;
  fl = [1 0; 0 1];
end

[Jdomain,JWdomain,qnode,qnode2,radius] = Jdomainf(tip,xyTip,enrich_node,4);

I1 = 0;
I2 = 0;
I  = [zeros(2,1)];
Iw  = [zeros(2,1)];
If  = [zeros(2,1)];
I_T  = [zeros(1,1)];
Iw_T  = [zeros(1,1)];
If_T  = [zeros(1,1)];

QT  = fl * [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];           % for the transformation to local coordinate
mu = E/(2.+ nu + nu);
kappa = 3-4*nu;    %Kolosov coeff, Plain strain

% ---------------------------
% Starting LOOP over ELEMENTS
%----------------------------
compt=0;
q=[];
done_nodes = [];

tstr = ['Tip',num2str(flag_end),': theta1 is ',num2str(t1),'   theta2 is ',num2str(t2),'   theta3  is ',num2str(t3),'\n'];
kstr = ['Tip',num2str(flag_end),': K1 is ',num2str(KI),'   K2 is ',num2str(KII),'   T is ',num2str(T),'  and theta is ',num2str(theta_inc),'\n'];
fprintf(output_file,kstr)
fprintf(output_file,tstr)



