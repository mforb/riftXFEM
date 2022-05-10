function [NP,NN,BP,BN] = nitcheNBmat(pt,pos,e,type_elem,enrich_node,xCrl,GVertex,crack_node,cont,pen)

%declare global variables here
global node element numnode numelem elemType
global incR xc yc phiN
global epsilon skip_branch

sctr = element(e,:);
nn   = length(sctr);
[N,dNdxi] = lagrange_basis(elemType,pt);  % element shape functions
J0 = node(sctr,:)'*dNdxi;                 % element Jacobian matrix
invJ0 = inv(J0);
dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
Gpt = N' * node(sctr,:);                  % GP in global coord, used

%Standard B matrix is computed always...
Nxfem = [ ] ;
Bxfem = [];
BxfemP = [];
BxfemN = [];
NxfemP = [];
NxfemN = [];
%loop on nodes, check is node is enriched........
if cont == 1
    Bfem = zeros(3,2*nn) ;
    Bfem(1,1:2:2*nn) = dNdx(:,1)' ;
    Bfem(2,2:2:2*nn) = dNdx(:,2)' ;
    Bfem(3,1:2:2*nn) = dNdx(:,2)' ;
    Bfem(3,2:2:2*nn) = dNdx(:,1)' ;
    Nfem = zeros(1,2*nn) ;
    Nfem(1,1:2:2*nn) = N' ;
    Nfem(1,2:2:2*nn) = N' ;
else
    Bfem = [] ;
    Nfem = [] ;
end

for in = 1:nn
    %Enriched by H(x) at global gauss point
    if ( enrich_node(sctr(in)) == 2)
        ref_elem = e;
        xCre = [xCrl(ref_elem,1) xCrl(ref_elem,2); xCrl(ref_elem,3) xCrl(ref_elem,4)];                 %each element has its crack!
        if ismember(sctr(in),crack_node) 
          Hi = sign(-1);
        else
          dist = signed_distance(xCre,node(sctr(in),:),0);
          Hi  = sign(dist);
        end
        % Bxfem at node "in"
        Hgp = 1;
        BI_enrP = [dNdx(in,1)*(Hgp - Hi) 0 ; 0 dNdx(in,2)*(Hgp - Hi) ;
            dNdx(in,2)*(Hgp - Hi) dNdx(in,1)*(Hgp - Hi)];
        N_enrP = [N(in)*(Hgp - Hi) , N(in)*(Hgp - Hi)] ;
        Hgp = -1;
        BI_enrN = [dNdx(in,1)*(Hgp - Hi) 0 ; 0 dNdx(in,2)*(Hgp - Hi) ;
            dNdx(in,2)*(Hgp - Hi) dNdx(in,1)*(Hgp - Hi)];
        N_enrN = [N(in)*(Hgp - Hi) , N(in)*(Hgp - Hi)] ;
        %else    %if the element is not cut by a crack, the enrichment is always 0 (NO LONGER TRUE)
            %BI_enr = [0 0 ; 0 0 ; 0 0];
        %end
        % Add to the total Bxfem
        BxfemP = [BxfemP BI_enrP];
        BxfemN = [BxfemN BI_enrN];
        NxfemP = [NxfemP N_enrP];
        NxfemN = [NxfemN N_enrN];
        clear BI_enrP; clear BI_enrN;
        clear N_enrP; clear N_enrN;
        
        % ------------ Enriched by asymptotic functions -----------------------------------
    elseif ( enrich_node(sctr(in)) == 1) % B(x) enriched node
        if type_elem(e,1) == 1   %looking for the "tip" element
            ref_elem = e;
        else    %trovo l'elemento/fessura a cui fa riferimento il nodo (SOLO 1 RIF AUTORIZZATO!!)
            [sctrn,xx] = find(element == sctr(in));
            [ele,xx] = find(type_elem(sctrn,:)==1);
            ref_elem = sctrn(ele);
            nR = find(enrich_node(sctr)==1);
            elem_blend = 1;
            Rpt = sum(N(nR));
            dRdx = sum(dNdx(nR,:),1);
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
        if ~elem_blend 
          aa = dNdx(in,1)*(Br(1)-BrI(1)) + N(in)*dBdx(1) ;
          bb = dNdx(in,2)*(Br(1)-BrI(1)) + N(in)*dBdy(1) ;
          B1_enr = [aa 0 ; 0 bb ; bb aa];
          
          aa = dNdx(in,1)*(Br(2)-BrI(2)) + N(in)*dBdx(2) ;
          bb = dNdx(in,2)*(Br(2)-BrI(2)) + N(in)*dBdy(2) ;
          B2_enr = [aa 0 ; 0 bb ; bb aa];
          
          aa = dNdx(in,1)*(Br(3)-BrI(3)) + N(in)*dBdx(3) ;
          bb = dNdx(in,2)*(Br(3)-BrI(3)) + N(in)*dBdy(3) ;
          B3_enr = [aa 0 ; 0 bb ; bb aa];
          
          aa = dNdx(in,1)*(Br(4)-BrI(4)) + N(in)*dBdx(4) ;
          bb = dNdx(in,2)*(Br(4)-BrI(4)) + N(in)*dBdy(4) ;
          B4_enr = [aa 0 ; 0 bb ; bb aa];
        else
          aa = Rpt*(dNdx(in,1)*(Br(1)-BrI(1)) + N(in)*dBdx(1)) + dRdx(1)*N(in)*(Br(1)-BrI(1)) ;
          bb = Rpt*(dNdx(in,2)*(Br(1)-BrI(1)) + N(in)*dBdy(1)) + dRdx(2)*N(in)*(Br(1)-BrI(1)) ;
          B1_enr = [aa 0 ; 0 bb ; bb aa];
          
          aa = Rpt*(dNdx(in,1)*(Br(2)-BrI(2)) + N(in)*dBdx(2)) + dRdx(1)*N(in)*(Br(2)-BrI(2)) ;
          bb = Rpt*(dNdx(in,2)*(Br(2)-BrI(2)) + N(in)*dBdy(2)) + dRdx(2)*N(in)*(Br(2)-BrI(2)) ;
          B2_enr = [aa 0 ; 0 bb ; bb aa];
          
          aa = Rpt*(dNdx(in,1)*(Br(3)-BrI(3)) + N(in)*dBdx(3)) + dRdx(1)*N(in)*(Br(3)-BrI(3)) ;
          bb = Rpt*(dNdx(in,2)*(Br(3)-BrI(3)) + N(in)*dBdy(3)) + dRdx(2)*N(in)*(Br(3)-BrI(3)) ;
          B3_enr = [aa 0 ; 0 bb ; bb aa];
          
          aa = Rpt*(dNdx(in,1)*(Br(4)-BrI(4)) + N(in)*dBdx(4)) + dRdx(1)*N(in)*(Br(4)-BrI(4)) ;
          bb = Rpt*(dNdx(in,2)*(Br(4)-BrI(4)) + N(in)*dBdy(4)) + dRdx(2)*N(in)*(Br(4)-BrI(4)) ;
          B4_enr = [aa 0 ; 0 bb ; bb aa];
        end
        if pen 
          aa = Rpt*N(in)*(Br(1))/sqrt(r);
          %aa = Rpt*N(in)*(Br(1))/r;
        else
          aa = Rpt*N(in)*(Br(1));
        end
        N1_enr = [aa 0 ; 0 aa];
        
        N_enr = [N1_enr];
        clear N1_enr;        
        Nxfem = [Nxfem N_enr];
        clear N_enr ;

        
        BI_enr = [B1_enr B2_enr B3_enr B4_enr];
        clear B1_enr; clear B2_enr; clear B3_enr; clear B4_enr;
        Bxfem = [Bxfem BI_enr];
        clear BI_enr ;
    end
end          % end of loop on nodes

BP = [ Bfem BxfemP];
BN = [ Bfem BxfemN ];
clear Bfem; clear BxfemP; clear BxfemN;
NP = [ Nfem NxfemP];
NN = [ Nfem NxfemN ];
clear Nfem; clear NxfemP; clear NxfemN;
