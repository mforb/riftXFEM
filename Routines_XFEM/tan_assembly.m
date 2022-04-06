function [sctrBfem,sctrBxfem] = tan_assembly(e,tan_element,pos)

global node orig_nn

sctr = tan_element(e,:);
nn   = length(sctr);
cnt = 0 ;

for k = 1 : nn
  cnt = cnt + 1 ;
  sctrBfem(2*k-1) = 2*sctr(k)-1 ;
  sctrBfem(2*k)   = 2*sctr(k)   ;
  sctrBxfem(2*cnt - 1) = 2 * pos(sctr(k)) - 1;
  sctrBxfem(2*cnt    ) = 2 * pos(sctr(k))    ;
end
