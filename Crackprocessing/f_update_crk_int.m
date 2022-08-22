function [ ckint_out, flag, crack_nodes ] = f_update_crk_int( ckint, intes, f1, f2, q1, q2, crack_nodes, sctr  )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Fri Nov 12 15:18:31 NZDT 2021

%this function makes sure that ckint is always a 2 point crack  in the same direction as the overall crack (q1 to q2)
% this function also returns a flag for there being overlapping intersects
% overlapping intersects mean the crack is going through a node
global node


flag = 0;
ins = intes(1,:);
if size(intes,1)>1
  for i = 2:size(intes,1)
    ft = 0;
    in = intes(i,:);
    for j = 1:size(ins,1)
      if points_same_2d(ins(j,:),in)
         flag = 1;
         ft = 1;
         for k = 1:length(sctr)
           if points_same_2d(in,node(sctr(k),:))
             crack_nodes = [ crack_nodes, sctr(k) ];
           end
         end
      end
    end
    if ft == 0
      ins = [ ins; in ];
    end
  end
end
       
crk = q2 - q1;

if isempty(ckint)
  if size(ins,1)==2
    ckint_out = f_align(ins(1,:),ins(2,:),crk);
  elseif (size(ins,1)==1) & not(f1 | f2)
    ckint_out = [q1, ins ];
  elseif f2
    ckint_out = [ ins , q2 ]; 
  elseif f1
    ckint_out = [ q1, ins ];
  end
end

if ~isempty(ckint)
  ckint_out = [ ckint(1,:),ins ];   % because we should always be going in the same direction n the rift
end
