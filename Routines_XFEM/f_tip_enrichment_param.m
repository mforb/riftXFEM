function [QT,tip,R,dRdx,Br,dBdx,dBdy] = f_tip_enrichment_param(e,Gpt,N,dNdx,sctr,xTip,xCrl,type_elem,enr_node,varargin)
global element 

f_op = 1;

in = find(enr_node(sctr) == 1,1);
if type_elem(e,1) == 1   %looking for the "tip" element
  ref_elem = e;
  R = 1;
  dRdx = [0, 0];
else    %trovo l'elemento/fessura a cui fa riferimento il nodo (SOLO 1 RIF AUTORIZZATO!!)
  [sctrn,xx] = find(element == sctr(in));
  [ele,xx] = find(type_elem(sctrn,:)==1);
  ref_elem = sctrn(ele);
  blend_elem = 1;
  nR = find(enr_node(sctr)==1);
  R = sum(N(nR));
  if ~isempty(dNdx)
    dRdx = sum(dNdx(nR,:),1);
  else
    dRdx = [0,0];
  end
end


if size(xTip,1)>1
  tip = xTip(ref_elem,:);
else
  tip = xTip;
end

if points_same_2d(xCrl(ref_elem,3:4),tip,1e-6)   
  xCrek  = [ xCrl(ref_elem,1:2); xCrl(ref_elem,3:4) ]; 
  seg   = xCrek(2,:) - xCrek(1,:);
else
  xCrek  = [ xCrl(ref_elem,3:4); xCrl(ref_elem,1:2) ]; 
  seg   = xCrek(1,:) - xCrek(2,:);
  f_op = -1;
end

alpha = atan2(seg(2),seg(1));
tip  = [xCrek(2,1) xCrek(2,2)];
QT    = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
xp    = QT*(Gpt-tip)';           % local coordinates

r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
theta = atan2(xp(2),xp(1));
if ( theta > pi | theta < -pi)
    disp (['something wrong with angle ',num2str(thet)]);
end

if abs(r) < 1e-8
  R = 0; 
end

if size(varargin,1)>0
  theta = pi*varargin{1};
  [Br,dBdx,dBdy] = branch_gp(r,theta,alpha);
else
  [Br,dBdx,dBdy] = branch_gp(r,theta,alpha);
end
