function [ Fc ] = f_crackforce( Fc, crack_lips, xCr, pos, enrich_node, split_elem, vertex_elem )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Mon Nov 29 11:58:33 NZDT 2021
global node element elemType

elems = union(split_elem,vertex_elem);

%NOT IMPLEMENTED YET

for kk = 1:size(xCr,2) %what's the crack?
  for ii=1:size(elems,1)
    iel = elems(ii) ;
    sctr=element(iel,:);
    p1 = crack_lips(ii,1:2,4,kk);
    p2 = crack_lips(ii,3:4,4,kk);
    p3 = crack_lips(ii,5:6,4,kk);
    d1 = crack_lips(ii,1:2,3,kk);
    d2 = crack_lips(ii,3:4,3,kk);
    d3 = crack_lips(ii,5:6,3,kk);
    dup1 = crack_lips(ii,1:2,1,kk);
    dup2 = crack_lips(ii,3:4,1,kk);
    dup3 = crack_lips(ii,5:6,1,kk);
    ddown1 = crack_lips(ii,1:2,2,kk);
    ddown2 = crack_lips(ii,3:4,2,kk);
    ddown3 = crack_lips(ii,5:6,2,kk);
    for k = 1:size(xCrk,2)
        sctrB = [sctrB assembly(iel,enrich_node(:,k),pos(:,k),k,crack_nodes)] ;
    end
    if ismember(iel,split_elem) 
          % here we are only going to be dealing in normal and tangential behaviour
          % we can work back to x-y after finding what the right forces are
       c_d1 = p1 + ddown1  
       c_p1 = p1 + dup1  
       c_d2 = p2 + ddown2
       c_p2 = p2 + dup2
       c_tan = c_d2 - c_d1;
       c_m   = c_d1 + 0.5*c_tan
       c_tan = c_tan/sqrt(c_tan(1)^2 +c_tan(2)^2);  
       c_norm = [ c_tan(2), -1*c_tan(1)];

    end
end

end
