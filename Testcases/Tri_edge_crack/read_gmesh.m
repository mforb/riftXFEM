%Square_mesh
iceshelfsquare

%[node,element,bcNodes,edgNodes]

node = msh.POS(:,1:2);
element = msh.TRIANGLES(:,1:3);

bs = find(node(:,2)==-1500);
[v,id] = unique(node(bs,1));
[v,ido] = sort(v,'ascend');

botNodes = bs(id(ido));
botEdges = [];
for i = 2:length(ido)
  botEdges = [ botEdges ; bs(id(ido(i-1))) bs(id(ido(i))) ]; 
end

rs = find(node(:,1)==1500);
[v,id] = unique(node(rs,2));
[v,ido] = sort(v,'ascend');

rightNodes = rs(id(ido));
rightEdges = [];
for i = 2:length(ido)
  rightEdges = [ rightEdges ; rs(id(ido(i-1))) rs(id(ido(i))) ]; 
end

ts = find(node(:,2)==1500);
[v,id] = unique(node(ts,1));
[v,ido] = sort(v,'descend');

topNodes = ts(id(ido));
topEdges = [];
for i = 2:length(ido)
  topEdges = [ topEdges ; ts(id(ido(i-1))) ts(id(ido(i))) ]; 
end

ls = find(node(:,1)==-1500);
[v,id] = unique(node(ls,2));
[v,ido] = sort(v,'descend');

leftNodes = ls(id(ido));
leftEdges = [];
for i = 2:length(ido)
 leftEdges = [ leftEdges ; ls(id(ido(i-1))) ls(id(ido(i))) ]; 
end


bcNodes = {botNodes rightNodes topNodes leftNodes}
edgNodes = {botEdges rightEdges topEdges leftEdges}
