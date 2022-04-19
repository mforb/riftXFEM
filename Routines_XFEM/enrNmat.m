function [Nxfem] = enrNmat(N,e,type_elem,enrich_node,xCrl,GVertex,cont)

%declare global variables here
global node element numnode numelem elemType
global incR xc yc phiN
global epsilon

sctr = element(e,:);
nn   = length(sctr);
Gpt = N' * node(sctr,:);                  % GP in global coord, used

%Standard B matrix is computed always...
Nxfem = [ ] ;
%loop on nodes, check is node is enriched........

for in = 1:nn
    %Enriched by H(x) at global gauss point
    if ( enrich_node(sctr(in)) == 2)
        ref_elem = e;
        N_enr = [N(in) 0 ; 0 N(in)] ;
        Nxfem = [Nxfem N_enr];
        clear N_enr ;
        
        % ------------ Enriched by asymptotic functions -----------------------------------
    elseif ( enrich_node(sctr(in)) == 1) % B(x) enriched node
        if type_elem(e,1) == 1   %looking for the "tip" element
            ref_elem = e;
        else    %trovo l'elemento/fessura a cui fa riferimento il nodo (SOLO 1 RIF AUTORIZZATO!!)
            [sctrn,xx] = find(element == sctr(in));
            [ele,xx] = find(type_elem(sctrn,:)==1);
            ref_elem = sctrn(ele);
        end
        % compute branch functions at Gauss point
        xCre  = [xCrl(ref_elem,1) xCrl(ref_elem,2); xCrl(ref_elem,3) xCrl(ref_elem,4)];
        seg   = xCre(2,:) - xCre(1,:);
        alpha = atan2(seg(2),seg(1));
        xTip  = [xCre(2,1) xCre(2,2)];
        QT    = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];
        xp    = QT*(Gpt-xTip)';           % local coordinates
        r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2));
        theta = atan2(xp(2),xp(1));
        
        if ( theta > pi | theta < -pi)
            disp (['something wrong with angle ',num2str(thet)]);
        end
        [Br,dBdx,dBdy] = branch_gp(r,theta,alpha);
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
        N1_enr = [aa 0 ; 0 aa];
        aa = N(in)*(Br(2)-BrI(2));
        N2_enr = [aa 0 ; 0 aa];
        aa = N(in)*(Br(3)-BrI(3));
        N3_enr = [aa 0 ; 0 aa];
        aa = N(in)*(Br(4)-BrI(4));
        N4_enr = [aa 0 ; 0 aa];
        
        
        
        N_enr = [N1_enr N2_enr N3_enr N4_enr];
        clear N1_enr; clear N2_enr; clear N3_enr; clear N4_enr;
        Nxfem = [Nxfem N_enr];
        clear N_enr ;
    elseif( enrich_node(sctr(in)) == 3)
        % compute the enrichment function and derivatives
        % at gauss point
        x = Gpt(1);
        y = Gpt(2);
        d = sqrt((x-xc)^2+(y-yc)^2);
        Fgp    = abs(d - incR) ;
        dFdx = sign(Fgp)*(x-xc)/d;
        dFdy = sign(Fgp)*(y-yc)/d;
        
        % compute the enrichment function at node
        x = node(sctr(in),1);
        y = node(sctr(in),2);
        d = sqrt((x-xc)^2+(y-yc)^2);
        FI = abs(d - incR);
        
        % Bxfem at node "in"
        aa = N(in)*(Fgp-FI) + N(in)*dFdx ;
        N_enr = [aa 0; 0 aa];
        
        % Add to the total Bxfem
        Nxfem = [Nxfem N_enr];
        clear N_enr ;
    end
    end          % end of loop on nodes
    % B matrix
    Nxfem;
end              % end of switch between enriched and non-enriched elements
%if e == 120
  %keyboard
%end
