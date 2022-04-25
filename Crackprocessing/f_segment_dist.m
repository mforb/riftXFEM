function [l,nv,mv,nnt,nmt,mmt] = f_segment_dist(seg)
x0 = seg(1) ; y0 = seg(2) ;
x1 = seg(3) ; y1 = seg(4) ;
l = sqrt((x1-x0)*(x1-x0) + (y1-y0)*(y1-y0)) ;
nv = [(y0-y1),(x1-x0)]./l;
mv = [(x1-x0),(y1 - y0)]./l;
nnt = nv'*nv;
nmt = mv'*nv;
mmt = mv'*mv;
