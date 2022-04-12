load ~/Data/RIS_Control_drag_april30mesh_iteration2_max210_MAC.mat
md=mechanicalproperties(md,md.results.StressbalanceSolution.Vx,md.results.StressbalanceSolution.Vy);

ISSM_xx = md.results.deviatoricstress.xx;
ISSM_yy = md.results.deviatoricstress.yy;
ISSM_xy = md.results.deviatoricstress.xy;
ISSM_h  = md.geometry.thickness;

grounded = md.mask.ocean_levelset;
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
element = elms(pk123,:);

TR = triangulation(element,node);
ext = freeBoundary(TR);
vx = md.initialization.vx;
vy = md.initialization.vy;
l = size(ext,1);
edges_front = [];
bc_fix = [];
for i = 1:l
  % check that the nodes are not in bc_front
  seg = ext(i,:); 
  t1 = ismember(seg(1),bc_front);
  t2 = ismember(seg(2),bc_front);
  if ~(t1&t2)
    % find normal to the segment
    x0 = node(seg(1),1) ; y0 = node(seg(1),2) ;
    x1 = node(seg(2),1) ; y1 = node(seg(2),2) ;
    vxl = mean([vx(seg(1)),vx(seg(2))]);
    vyl = mean([vy(seg(1)),vy(seg(2))]);
    l = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0)) ;
    ns = [(y0-y1),(x1-x0)]./l;
    vn = ns*[vxl;vyl];
    if abs(vn) < 10 
      bc_fix = [bc_fix; seg(1) ; seg(2) ];
    end
  else
    edges_front = [edges_front; seg ];
  end
end
bc_fix = unique(bc_fix);
    

mv = intersect(pk1,pk2);
mid = intersect(mv,pk3);
edg_elm = setdiff(pk123,mid);
edg_elm_nodes = unique(elms(edg_elm,:));
size(edg_elm_nodes)


ISSM_xx = ISSM_xx(pk123);
ISSM_yy = ISSM_yy(pk123);
ISSM_xy = ISSM_xy(pk123);
ISSM_H = [];
for i = 1:length(element)
  sctr = element(i,:);
  ISSM_H(i) = mean(ISSM_h(sctr));
end

save('/home/antarctica/Data/import_issm_holly1.mat','element','node','bc_fix','bc_front','edges_front','ISSM_xx','ISSM_yy','ISSM_xy','ISSM_H');
