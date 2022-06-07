function [N,Gpt] = xfemNmat(pt,e,type_elem,enr_node,xCrl,GVertex,xTip,crack_node,cont,varargin)

%declare global variables here
global node element numnode numelem elemType
global incR xc yc phiN
global epsilon

sctr = element(e,:);
nn   = length(sctr);
[N,dNdxi] = lagrange_basis(elemType,pt);  % element shape functions
Gpt = N' * node(sctr,:);                  % GP in global coord, used
xCre = xCrl(e,:);                 %each element has its crack!
if size(varargin,1)>0
  Hgp = varargin{1};
else
  dist = signed_distance(xCre,Gpt,0);
  Hgp  = sign(dist);
end

%Standard B matrix is computed always...
if cont == 1
    Nfem = zeros(2,2*nn) ;
    Nfem(1,1:2:2*nn) = N' ;
    Nfem(2,2:2:2*nn) = N' ;
else
    Bfem = [] ;
end
blend_elem = 0;
if any(enr_node(sctr) == 1)
  in = find(enr_node(sctr) == 1,1);
  if type_elem(e,1) == 1   %looking for the "tip" element
    ref_elem = e;
    Rpt = 1;
  else    %trovo l'elemento/fessura a cui fa riferimento il nodo (SOLO 1 RIF AUTORIZZATO!!)
    [sctrn,xx] = find(element == sctr(in));
    [ele,xx] = find(type_elem(sctrn,:)==1);
    ref_elem = sctrn(ele);
    blend_elem = 1;
    nR = find(enr_node(sctr)==1);
    Rpt = sum(N(nR));
  end
  if points_same_2d(xCrl(ref_elem,3:4),xTip(ref_elem,:),1e-6)   
    xCrek  = [ xCrl(ref_elem,1:2); xCrl(ref_elem,3:4) ]; 
  else
    xCrek  = [ xCrl(ref_elem,3:4); xCrl(ref_elem,1:2) ]; 
  end
  seg   = xCrek(2,:) - xCrek(1,:);
  alpha = atan2(seg(2),seg(1));
  xTip  = [xCrek(2,1) xCrek(2,2)];
  QT    = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
  xp    = QT*(Gpt-xTip)';           % local coordinates
  if abs(xp) < 1e-6
    Rpt = 0; % the point is (in theory) the tip, so the enrichments should all be zero
  end
  r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
  theta = atan2(xp(2),xp(1));
  if ( theta > pi | theta < -pi)
      disp (['something wrong with angle ',num2str(thet)]);
  end
  if size(varargin,1)>0
    theta = pi*varargin{1};
    [Br,dBdx,dBdy] = branch_gp(r,theta,alpha);
  else
    [Br,dBdx,dBdy] = branch_gp(r,theta,alpha);
  end
end

iwant = [ ] ;
elem_blend = 0;
%switch between non-enriched and enriched elements
if( all(enr_node(sctr) == 0) ) 
  N = Nfem ;
else
  Nxfem = [ ] ;
  %loop on nodes, check is node is enriched........
   
  for in = 1:nn
    %Enriched by H(x) at global gauss point
    if ( enr_node(sctr(in)) == 2)
      if (type_elem(e,cont) == 3)
        if isempty(varargin)
          distV = signed_distance(xCre,GVertex(e,:),0);
          HV = sign(distV);
          if HV * Hgp <= 0 %below or above the line relating the crack intersections?
              Hgp = Hgp;
          elseif abs(distV) < epsilon
              Hgp = Hgp;
          else
              vv = [xCre(1,:); xCre(2,:); GVertex(e,:)];
              flag = inhull(Gpt,vv); % is the point in the triangle formed by the crack ends and the vertex?
              if flag == 1     % yes
                  Hgp =- Hgp;
              else             % no
                  Hgp =  Hgp;
              end
          end
        end
        dist = signed_distance(xCre,node(sctr(in),:),0);
        Hi   = sign(dist);
        if ismember(sctr(in),crack_node) 
          %Hi = sign(-1);
          %Hi = -1*Hgp;
          Hi = 0;
        end
      else%
        % Enrichment function, H(x) at node "in"
        dist = signed_distance(xCre,node(sctr(in),:),0);
        Hi  = sign(dist);
        if ismember(sctr(in),crack_node) 
          %Hi = sign(-1);
          %Hi = -1*Hgp;
          Hi = 0;
        end
      end % is vertex
      % Bxfem at node "in"
      NI_enr = [N(in)*(Hgp - Hi) 0 ; 0 N(in)*(Hgp - Hi) ];
      % Add to the total Nxfem
      Nxfem = [Nxfem NI_enr];
      clear NI_enr ;
        
    % ------------ Enriched by asymptotic functions -----------------------------------
    elseif ( enr_node(sctr(in)) == 1) % B(x) enriched node

      
      % compute branch functions at node "in"
      xp    = QT*(node(sctr(in),:)-xTip)';
      r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
      theta = atan2(xp(2),xp(1));
      
      if ( theta > pi | theta < -pi)
          disp (['something wrong with angle ',num2str(thet)]);
      end
      [BrI] = branch_node(r,theta);
      
      % composants of Benr matrix
      aa = N(in)*(Br(1)-BrI(1));
      N1_enr = [aa 0 ; 0 aa ];
      
      aa = N(in)*(Br(2)-BrI(2));
      N2_enr = [aa 0 ; 0 aa];
      
      aa = N(in)*(Br(3)-BrI(3));
      N3_enr = [aa 0 ; 0 aa];
      
      aa = N(in)*(Br(4)-BrI(4)) ;
      N4_enr = [aa 0 ; 0 aa];
      
      NI_enr = Rpt*[N1_enr N2_enr N3_enr N4_enr];
      clear N1_enr; clear N2_enr; clear N3_enr; clear N4_enr;
      Nxfem = [Nxfem NI_enr];
      clear NI_enr ;
    end
  end          % end of loop on nodes
  % B matrix
  N = [ Nfem Nxfem ];
  clear Nfem; clear Nxfem;
end              % end of switch between enriched and non-enriched elements
