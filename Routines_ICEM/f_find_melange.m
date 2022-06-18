function [flag1,width,nodeTanfix] = f_find_melange(e,xCrk,nodeTanfix);
global node element 
flag1 = 0
for kj = 1:size(xCrk.coor,1)-1       %loop over the elements of the fracture
  if xCrk.melange(kj)
    q1 = xCrk(kk).coor(kj,:); 
    q2 = xCrk(kk).coor(kj+1,:);
    [f1,f2,f3,cn,~,d_r] = crack_interact_element([q1,q2],e,[]);
    if f1
      flag1 = 1;
      if flag3  
        width = xCrk.width(kj+1);
        break
      else
        width = [mel_elems; e xCrk.width(kj)*d_r(1)+xCrk.width(kj+1)*d_r(2) ];
      end
    elseif ~isempty(cn)
      flag1 = 1;
      p = mean(cn,1);
      d1 = norm(q1-p);
      d2 = norm(q2-p);
      width = (xCrk.width(kj)*d1+xCrk.width(kj+1)*d2)/2*(d1+d2);
      if nargin > 2 
        sctr = element(e,:);
        for i=1:length(sctr)
          nod = sctr(i);
          nodeTanfix(nodeTanfix==nod) = [];
        end
      end
      break
    end
end
