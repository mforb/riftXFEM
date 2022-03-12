
load ~/Data/Holly/RIS_Control_drag_april30mesh_iteration2_max210_MAC.mat

grounded = md.mask.groundedice_levelset;
ice      = md.mask.ice_levelset;

elms = md.mesh.elements;
node  = [md.mesh.x,md.mesh.y];

bc_front =find(ice==0);
n_keep  = find(grounded<0);

v1 = intersect(elms(:,1),n_keep);
v2 = intersect(elms(:,2),n_keep);
v3 = intersect(elms(:,3),n_keep);

pk1 = find(ismember(elms(:,1),v1));
pk2 = find(ismember(elms(:,1),v2));
pk3 = find(ismember(elms(:,1),v3));

pk12 = union(pk1,pk2);
pk123 = union(pk12,pk3);
elements = elms(pk123,:);

mv = intersect(pk1,pk2);
mid = intersect(mv,pk3);
edg_elm = setdiff(pk123,mid);
edg_elm_nodes = unique(elms(edg_elm,:));
size(edg_elm_nodes)
bc_fix = setdiff(edg_elm_nodes,n_keep);

element = elms(pk123,:);

