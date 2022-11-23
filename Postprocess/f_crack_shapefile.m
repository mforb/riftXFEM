function [] = f_crack_shapefile(xCr,results_path,sname)

sr = struct();
sr.Geometry = 'Line';
sr.BoundingBox = [min(xCr.coor(:,1)), min(xCr.coor(:,2)); max(xCr.coor(:,1)), max(xCr.coor(:,2))];
sr.X = xCr.coor(:,1)';
sr.Y = xCr.coor(:,2)';
sr.id = NaN;
shapewrite(sr,[results_path,'/',sname]);

if isfield(xCr,'melange')
  sp = polyshape();
  inds = find(xCr.melange)
  for i = 1:length(inds)
    id = inds(i);
    wid = [xCr.width(id),xCr.width(id+1)]
    sv = xCr.coor(id+1,:)-xCr.coor(id,:)
    nv = [sv(2),-sv(1)]/(norm(sv));
    c1 = xCr.coor(id,:)+0.5*wid(1)*nv;
    c2 = xCr.coor(id,:)-0.5*wid(1)*nv;
    c3 = xCr.coor(id+1,:)-0.5*wid(2)*nv;
    c4 = xCr.coor(id+1,:)+0.5*wid(2)*nv;
    cs = [c1;c2;c3;c4];
    ptemp = polyshape(cs(:,1),cs(:,2));
    if i == 1
      sp = ptemp;
    else
      sp = union(sp,ptemp);
    end
  end
  %keyboard

  sr = struct();
  sr.Geometry = 'Polygon';
  sr.BoundingBox = [min(sp.Vertices(:,1)), min(sp.Vertices(:,2)); max(sp.Vertices(:,1)), max(sp.Vertices(:,2))];
  sr.X = sp.Vertices(:,1)';
  sr.Y = sp.Vertices(:,2)';
  sr.id = NaN;
  shapewrite(sr,[results_path,'/',sname,'_mp']);
  %sp.Geometry = 'PolyLine'
  %sp.id = NaN
end





