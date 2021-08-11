function [enrdomain,radius] = ENRdomainf(X0,X1)

global node element

numnode = size(node,1);
% -------------------------------------
radius = sqrt((X1(1)-X0(1))^2+(X1(2)-X0(2))^2);
center = X0;

r=[];
% Distance from the center of tip element
for i = 1 : numnode
    sctr = node(i,:);
    rho  = sqrt((sctr(1)-center(1))^2+(sctr(2)-center(2))^2) ;
    r    = [r,rho];
end

test = r-radius ;
test = test(element)';
test = min(test);
enrdomain = find(test<=0)';



