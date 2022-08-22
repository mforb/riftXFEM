function [Nxfem] = enrNmat(N,e,type_elem,enr_node,xCrl,GVertex,xTip,cont,pen)

%declare global variables here
global node element numnode numelem elemType
global incR xc yc phiN
global epsilon skip_branch

sctr = element(e,:);
nn   = length(sctr);
Gpt = N' * node(sctr,:);                  % GP in global coord, used

%Standard B matrix is computed always...
Nxfem = [ ] ;
%loop on nodes, check is node is enriched........
if any(enr_node(sctr) == 1)
    [QT,tip,Rpt,~,Br,~,~] = f_tip_enrichment_param(e,Gpt,N,[],sctr,xTip,xCrl,type_elem,enr_node,1); % this is always the "positive" sideas defined by the crack orientation 
end

for in = 1:nn
    %Enriched by H(x) at global gauss point
    if ( enr_node(sctr(in)) == 2)
        ref_elem = e;
        N_enr = [N(in) 0 ; 0 N(in)] ;
        Nxfem = [Nxfem N_enr];
        clear N_enr ;
        
        % ------------ Enriched by asymptotic functions -----------------------------------
    elseif ( enr_node(sctr(in)) == 1) % B(x) enriched node
       if ~skip_branch
          % compute branch functions at node "in"
          if pen 
            aa = Rpt*N(in); % (Br(1))/sqrt(r) should equal 1 in this case
            %aa = Rpt*N(in)*(Br(1))/r;
          else
            aa = Rpt*N(in)*(Br(1));
          end
          N1_enr = [aa 0 ; 0 aa];
          
          N_enr = [N1_enr];
          clear N1_enr;        
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
