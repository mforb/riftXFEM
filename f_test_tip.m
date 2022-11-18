function [ tt ] = f_test_tip( tip, delta_inc)
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
global node element elemType

% This function determins if a tip is within an element with sufficient margin 
% first thing is just to check if any noodes are within epsilon 
tt = 0
x = tip(1);
y = tip(2);
p_buff1 = polybuffer(tip, 'points',delta_inc);
inds = find( isinterior(p_buff1, nodes(:,1), nodes(:,2)) ) ;
nds = node(inds,:)

p_buff = polybuffer(tip, 'points',epsilon);
if any( isinterior(p_buff,nds(:,1),nds(:,2)) )
  tt = 1 
  return
end

% find closest point
ndst_x = nds(:,1) - tip(1); 
ndst_y = nds(:,2) - tip(2); 
ndst = ndst_x.*ndst_x + ndst_y.*ndst_y
in = find(ndst == min(ndst));
[els,~] = find(element==nds(in));
for i = 1:length(els)
  e = els(i);
  sctr=element(e,:);
  sctrl = [sctr sctr(1,1)];      
  for iedge=1:size(sctr,2)             %loop over the edges of elements
    nnode1=sctrl(iedge);
    nnode2=sctrl(iedge+1);
    p1 = [node(nnode1,:)];
    p2 = [node(nnode2,:)];
    y0 = p1(2); x0 = p1(1);
    y1 = p2(2); x1 = p2(1);
    l = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0)) ;

    psi1 = (x-x1)*(x1-x0)+(y-y1)*(y1-y0);
    psi = psi/l;

    psi2 = (x-x0)*(x0-x1)+(y-y0)*(y0-y1);
    psi = psi/l;
    
    phi = (y0-y1)*(x-x0) + (x1-x0)*(y-y0);
    phi = phi/l;
    if phi < epsilon 
    if psi1 < 0 & psi2 < 0 
      tt = 1
      return
    end
  end
end
