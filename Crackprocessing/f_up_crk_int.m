function [ ckint_out ] = f_up_crk_int( ckint, intes, f1, f2, q1, q2, el  )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Fri Nov 12 15:18:31 NZDT 2021

%this function makes sure that ckint is always a 2 point crack  in the same direction as the overall crack (q1 to q2)
global node epsilon

if epsilon > 1e-6
  in_epsi = 1e-6;
else
  in_epsi = epsilon/10; 
end

ins = intes(1,:);
if size(intes,1)>1
for i = 2:size(intes,1)
    skip = 0;
    in = intes(i,:);
    for j = 1:size(ins,1)
      if points_same_2d(ins(j,:),in,in_epsi) % the stricter confidence is to protect long skinny elements, or vertex edge proximity from stopping the calculation
         skip = 1;
      end
    end
    if ~skip 
      ins = [ ins; in ];
    end
end
end

if size(ins,1)>2
  error(['There are more then 2 intersects in element ',num2str(el)])
end
       
crk = q2 - q1;

if isempty(ckint)
  if size(ins,1)==2
    ckint_out = f_align(ins(1,:),ins(2,:),crk);
  elseif f2
    ckint_out = [ ins , q2 ]; 
  elseif f1
    ckint_out = [ q1, ins ];
  else
    msg = ['Crk_int upgrade error, element ',num2str(el)]
    error(msg)
  end
else
  if f1
    ckint_out = [ ckint(1:2),ins(1,:) ];   % because we should always be going in the same direction n the rik
   % there should be only one ins at this point
  else
    msg = ['Crk_int upgrade error, element ',num2str(el)]
    error(msg)
  end

end
