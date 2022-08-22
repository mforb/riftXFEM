function [ckint,flag] = update_crk_int(ckint,intes,q1,q2)
flag = 0;
if ~isempty(ckint) 
  N = length(ckint)/2
  for n = 1:N 
    if points_same_2d([ckint(n),ckint(n+1)],pt)  
      flag = 1; 
    end
  end
end
if ~flag
  ckint = [ckint,pt];
end
end
