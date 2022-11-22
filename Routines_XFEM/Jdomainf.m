function [Jdomain,JWdomain,qnode,qnode2,radius] = Jdomainf(tip_elem,xTip,enrich_node,fac)

global node element typeProblem

numnode = size(node,1);

% calculation of the area of the tip element
% x = node(element(tip_elem,:),:);
% % Area = sum of areas of each sub-triangle
% x0 = x(1,1);
% y0 = x(1,2);
% 
% x1 = x(2,1);
% y1 = x(2,2);
% 
% x2 = x(3,1);
% y2 = x(3,2);
% 
% x3 = x(4,1);
% y3 = x(4,2);
% 
% A1 = 0.5 * ((x0-x2)*(y1-y2) - (x1-x2)*(y0-y2)) ;
% A2 = 0.5 * ((x0-x3)*(y2-y3) - (x2-x3)*(y0-y3)) ;
% area = A1 + A2 
sctr = element(tip_elem,:);
area = polyarea(node(sctr,1),node(sctr,2));


% J radius = fac * sqrt(area);
%if strcmp(typeProblem,'ISSM')
  %fac = 4;
%else
  %fac = 3;
%end
%fac = 6;
radius = fac * sqrt(area);
center = xTip;

r=[];
% Distance from the center of tip element
for i = 1 : numnode
    sctr = node(i,:);
    rho  = sqrt((sctr(1)-center(1))^2+(sctr(2)-center(2))^2);
    r    = [r,rho];
end
test = r-radius;
test = test(element)';
testW = min(test);
test = max(test).*min(test);
Jdomain = find(test<=0);
JWdomain = find(testW<0);
test1 = r-radius;
test2 = test1(element(JWdomain,:))';
test1 = test1(element(Jdomain,:))';
test1 = (test1<=0);
test2 = test2<=0;
qnode = test1';
qnode2 = test2';

% checking if there is some B enriched elements in the Jdomain
for i=1:size(Jdomain,2)
    sctr = element (Jdomain(i),:);
    tip_enr = find(enrich_node(sctr,:) == 1);
    if size(tip_enr,1) > 0 
        disp('!!!WARNING!!! the Jdomain contain tip enriched elements => you should use a finer mesh')
        return % in order to not display the error message several times 
    end
end

