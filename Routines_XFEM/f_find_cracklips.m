function [ crack_lips, Flag_pen ] = f_find_cracklips( u, xCr, kk, enr_dom, type_elem, xCrl,xTip,xVertex,enr_node,crack_node,pos,split_elem, vertex_elem, tip_elem)
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Thu Nov 25 11:49:40 NZDT 2021
global node element elemType
global E nu C sigmato
global Jint iMethod
global epsilon

if strcmp(elemType,'Q4')
  corner = [1 2 3 4 1] ;
  nnode = [-1 -1;1 -1;1 1;-1 1] ;
elseif strcmp(elemType, 'T3')
  corner = [1 2 3 1] ;
  nnode = [0 0;1 0;0 1] ;
end

elems = union(split_elem,vertex_elem);
elems = union(elems,tip_elem);

crack_lips = zeros( size(elems,1),6,4,size(xCr,2));
Flag_pen = 0;

for ii=1:length(elems)
  p = [];
  iel = elems(ii) ;
  sctr=element(iel,:);
  nn = length(sctr);
  vv = node(sctr,:);
  [phi] = dista(iel,xCrl) ;
  % First we find out what the points are that we are going to evalaute the crack lips at
  if ismember(iel, tip_elem)
    tip = xTip(iel,:);
    ntip = f_naturalpoint(tip,vv,20,1e-6);
    psi = f_dista2(iel,xCrl,tip);
    inter = intersect(sctr,crack_node);
    if any(inter);
      int = find(sctr==inter);
      p = [nnodes(int,:) ; ntip ];
    else
      [cutEdge,nnodes] = f_edgedetect(nnode, corner,  phi, psi) ;
      p = [nnodes(end,:) ; ntip ];
    end
  else 
    [cutEdge, nnodes] = f_edgedetect(nnode, corner,  phi) ;
    nEdge = length(cutEdge);
    if nEdge==1 % then one of the nodes must be a crack_node
      crack_n = intersect(crack_node,sctr);
      if isempty(crack_n)
        error("only one edge but no crack node in elem ,",num2str(iel),", when evaluating crack lips")
      end
      crack_c = find(sctr==crack_n);
      nnodes = [ nnodes; nnodes(crack_c,:) ] ;
    end

    p = [nnodes(end-1,:);nnodes(end,:)];

    if ismember(iel,vertex_elem)
      tip = xVertex(iel,:);
      ntip = f_naturalpoint(tip,vv,20,1e-6);
      p = [p(1,:) ; ntip ; p(2,:) ];
    end
  end

  % now we are going to evaluate the displacement of each point for the positive, negative and midpoint of the crack
 %[A,BrI,QT,Tip,alpha] = f_enrich_assembly(iel,pos,type_elem,xCrl,enr_node);
 AB = assembly(iel,enr_node,pos,kk,crack_node);
 uAB = u(AB); 

 for gp = 1:size(p,1) 
   gpt = p(gp,:);
   c_inds = [2*gp - 1, 2*gp ];
   Hgp = 1;
   [Nmat_top,Gpt] = xfemNmat(gpt,iel,type_elem,enr_node(:,kk),xCrl,xVertex,xTip,crack_node,1,Hgp);
   Hgp = -1;
   [Nmat_bot,~] = xfemNmat(gpt,iel,type_elem,enr_node(:,kk),xCrl,xVertex,xTip,crack_node,1,Hgp);
   Nmid = Nmat_top; 
   Nmid(:,(2*nn+1):end) = 0;
   Ntop = Nmat_top;
   Ntop(:,1:2*nn) = 0;
   Nbot = Nmat_bot;
   Nbot(:,1:2*nn) = 0;
   crack_lips(ii,c_inds,1,kk) = Gpt ;
   crack_lips(ii,c_inds,2,kk) = Ntop*uAB;
   crack_lips(ii,c_inds,3,kk) = Nbot*uAB;
   crack_lips(ii,c_inds,4,kk) = Nmid*uAB;
   %if ismember(iel, tip_elem)
     %keyboard
   %end
   % check that there is no interpenetration
   [~,nv,~,~,~,~] = f_segment_dist(xCrl(iel,:));
   gn = nv*(Ntop*uAB - Nbot*uAB);
   if gn < 0 % crack tips can be a problem for this 
     Flag_pen = 1;
   end
 end
end
