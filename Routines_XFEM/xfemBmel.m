function [B,dJO] = xfemBmel(pt,e,type_elem,enrich_node,xCrl,GVertex,crack_node,cont,QT,mT)


%declare global variables here
global node element numnode numelem elemType
global incR xc yc phiN
global epsilon

vN = QT(:,2)';

sctr = element(e,:);
nn   = length(sctr);
[N,dNdxi] = lagrange_basis(elemType,pt);  % element shape functions
J0 = node(sctr,:)'*dNdxi;                 % element Jacobian matrix
nodeN = node(sctr,:)*vN';            % nodes in coordinate system perp
rN = max(nodeN)-min(nodeN);                        % distance between nodes and this is what we want to use to stretch/shrink the elements transform
invJ0 = inv(J0);
dNdx  = dNdxi*invJ0;                      % derivatives of N w.r.t XY
% now we want to rotate using QT
dNdxN = QT*dNdx';
dNdxN(2,:) = dNdxN(2,:)*rN/mT;
dNdx = ( QT'*dNdxN )';

dJO = det(J0)*rN/mT;

Gpt = N' * node(sctr,:);                  % GP in global coord, used

%%Standard B matrix is computed always...
%if cont == 1
    %Bfem = zeros(3,2*nn) ;
    %Bfem(1,1:2:2*nn) = dNdx(:,1)' ;
    %Bfem(2,2:2:2*nn) = dNdx(:,2)' ;
    %Bfem(3,1:2:2*nn) = dNdx(:,2)' ;
    %Bfem(3,2:2:2*nn) = dNdx(:,1)' ;
%else
    %Bfem = [] ;
%end

Bfem = [] ;

%iwant = [ ] ;
%switch between non-enriched and enriched elements
if( any(enrich_node(sctr)) == 0 ) & isempty(intersect(crack_node,sctr)) 
    %if (type_elem(e,cont) == 4 )
      %ref_elem = e;
      %xCre = [xCrl(ref_elem,1) xCrl(ref_elem,2); xCrl(ref_elem,3) xCrl(ref_elem,4)];                 %each element has its crack!

    %else
    warning('Enrichment 0 in melange set-up')
else
    Bxfem = [ ] ;
    %loop on nodes, check is node is enriched........
    
    for in = 1:nn
        %Enriched by H(x) at global gauss point
        if ( enrich_node(sctr(in)) == 2)
            ref_elem = e;
            xCre = [xCrl(ref_elem,1) xCrl(ref_elem,2); xCrl(ref_elem,3) xCrl(ref_elem,4)];                 %each element has its crack!
            if (type_elem(e,cont) == 3)
                distV = signed_distance(xCre,GVertex(ref_elem,:),0);
                HV = sign(distV);
                dist = signed_distance(xCre,Gpt,0);
                Hgp  = sign(dist);
                if HV * Hgp <= 0 %below or above the line relating the crack intersections?
                    Hgp = Hgp;
                elseif abs(distV) < epsilon
                    Hgp = Hgp;
                else
                    vv = [xCre(1,:); xCre(2,:); GVertex(ref_elem,:)];
                    flag = inhull(Gpt,vv); % is the point in the triangle formed by the crack ends and the vertex?
                    if flag == 1     % yes
                        Hgp =- Hgp;
                    else             % no
                        Hgp =  Hgp;
                    end
                end
                dist = signed_distance(xCre,node(sctr(in),:),0); % nodes are always outside the triangle..easy!
                Hi  = sign(dist);
                
            else%
                % Enrichment function, H(x) at global Gauss point
                dist = signed_distance(xCre,Gpt,16);
                Hgp  = sign(dist);
                % Enrichment function, H(x) at node "in"
                if ismember(sctr(in),crack_node) 
                  Hi = sign(-1);
                else
                  dist = signed_distance(xCre,node(sctr(in),:),0);
                  Hi  = sign(dist);
                end
            end % is vertex
            % Bxfem at node "in"
            %BI_enr = [dNdx(in,1)*Hi 0 ; 0 dNdx(in,2)*Hi ; dNdx(in,2)*Hi dNdx(in,1)*Hi];
            BI_enr = [dNdx(in,1)*(Hgp - Hi) 0 ; 0 dNdx(in,2)*(Hgp - Hi) ;
                dNdx(in,2)*(Hgp - Hi) dNdx(in,1)*(Hgp - Hi)];
            %BI_enr = [dNdx(in,1) 0 ; 0 dNdx(in,2) ; dNdx(in,2) dNdx(in,1)];
            %else    %if the element is not cut by a crack, the enrichment is always 0 (NO LONGER TRUE)
                %BI_enr = [0 0 ; 0 0 ; 0 0];
            %end
            % Add to the total Bxfem
            Bxfem = [Bxfem BI_enr];
            clear BI_enr ;
            
            % ------------ Enriched by asymptotic functions -----------------------------------
        elseif ( enrich_node(sctr(in)) == 1) % B(x) enriched node
          warning('Enrichment 1 in melange set-up')
        end
    end          % end of loop on nodes
    % B matrix
    B = [ Bfem Bxfem ];
    clear Bfem; clear Bxfem;
end              % end of switch between enriched and non-enriched elements
%if e == 120
  %keyboard
%end
