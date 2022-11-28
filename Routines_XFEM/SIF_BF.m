function [Knum,theta_inc] = SIF_BF(C,flag_end,tip,elem_crk,xCr,type_elem,enrich_node,crack_nodes,xVertex,pos,u,F,kk,alpha,...
    tip_elem,split_elem,vertex_elem,corner_elem,tan_elem,elem_force)
    % main difference here is that F is also passed into the function

global node element elemType E nu Cm1
global iMethod iParam
global incR xc yc phiN
global lambda1 lambda2 nu1 nu2
global elemType typeMesh typeProblem typeCrack stressState
global wall_int output_file wall_force crack_load fixwf
global quick_freeze melange

plotNode = 'NO' ;

if isempty(wall_force)
  wall_force = 1;
end
if ~isempty(crack_load) & crack_load~=0
  wall_force = 1;
end
if melange
  wall_force = 1;
  elems = union(split_elem,vertex_elem);
  elemst = tan_elem;

  mel_elems = [];
  tot_phi = 0;
  for kk = 1:size(xCr,2)
    for i=1:length(elems)                     %loop on elems (=elements selected for enrichment)
      iel = elems(i);
      [flag1,width,phiR,~] = f_find_melange(iel,xCr(kk),[]);
      if flag1
        tot_phi = tot_phi + phiR;
        mel_elems = [mel_elems; kk, iel, width];
      end
    end
  end
  mE = tot_phi/size(mel_elems,1);
end

if ~isempty(fixwf) & fixwf == 1 
  wall_force = 0  
end

if strcmp(elemType,'Q4') 
  intType = 'GAUSS' ;
  corner = [1 2 3 4 1] ;
  nnode = [-1 -1;1 -1;1 1;-1 1] ;
else
  intType = 'TRIANGULAR';
  corner = [1 2 3 1] ;
  nnode = [0 0;1 0;0 1] ;
end
% for debug
%TR = triangulation(element,node);
%figure(12)
%triplot(TR)
%hold on

% Compute the Stress Intensity Factors
% Using the Interaction integral method

% Steps :
% 1- detection of the elements on which we integrate
% 2- loop over these elements
% 3- loop over Gauss points
% 4- computation of stress, strain... in local coordinates !!!   ATTENTION
% 5- computation of the auxilliary fields: AuxStress and AuxEps and AuxGradDisp
% 6- computation of I1 and I2

% Determine J domain and weight function
if flag_end == 1
  xyTip = [elem_crk(tip,1) elem_crk(tip,2)] ;
  fl = [1 0; 0 -1];
elseif flag_end == 2
  xyTip = [elem_crk(tip,3) elem_crk(tip,4)] ;
  fl = [1 0; 0 1];
end

[Jdomain,JWdomain,qnode,qnode2,radius] = Jdomainf(tip,xyTip,enrich_node,5);

I1 = 0;
I2 = 0;
I  = [zeros(2,1)];
Iw  = [zeros(2,1)];
If  = [zeros(2,1)];

QT  = fl * [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];           % for the transformation to local coordinate
mu = E/(2.+ nu + nu);
kappa = 3-4*nu;    %Kolosov coeff, Plain strain

% ---------------------------
% Starting LOOP over ELEMENTS
%----------------------------
compt=0;
q=[];
done_nodes = [];

for iel = 1 : size(JWdomain,2)
    e = JWdomain(iel) ; % current element
    sctr = element(e,:);
    nn   = length(sctr);
    normal_order = 6;          % max = 8
    tip_order    = 20;         % is not used -> if >20 allow a error flag (bug) if the J integrale is performed on blending's
    split_order  = 7;
    vertex_order = 7;
    qf     = qnode2(iel,:);
    mel_bool = 0;

    if ismember(e,Jdomain)
      if( ismember(e,split_elem) && any(ismember(enrich_node(sctr),3)) ) %why?
          order = 3 ;
          phis = phiN(sctr) ;
          if strcmp(elemType,'Q4') 
            [W,Q] = discontQ4quad(order,phis,nnode) ;
          else
            [W,Q] = discontT3(order,phis,nnode) ;
          end
      else
          %choose Gauss quadrature rules for elements
          [W,Q] = gauss_rule(e,enrich_node,elem_crk,...
              xyTip,xVertex,tip_elem,split_elem,vertex_elem,corner_elem,xCr) ;
      end
      
      %     compt=compt+1;
      %     plot_GP_new(elemType,q,Q,xCr,W,e,compt,enrich_node,[])
      
      % ----------------------------
      % start loop over Gauss points
      % -----------------------------
      for q = 1:size(W,1)
          pt = Q(q,:);
          wt = W(q);
          [N,dNdxi] = lagrange_basis(elemType,pt);
          J0    = node(sctr,:)'*dNdxi;
          invJ0 = inv(J0);
          dNdx  = dNdxi*invJ0;
          gradq = qf*dNdx;
          Gpt = N' * node(sctr,:);     % GP in global coord
          B = [];
          % +++++++++++++++++++++++++
          % Gradient of displacement
          % +++++++++++++++++++++++++
          % need to compute u,x u,y v,x v,y, stored in matrix H
          for k = 1:size(xCr,2)
              B = [B xfemBmat(pt,e,type_elem,enrich_node(:,k),elem_crk,xVertex,xyTip,crack_nodes,k)];
          end

          leB = size(B,2);
          
          % nodal displacement of current element
          % taken from the total nodal parameters u
          U = [];
          for k = 1:size(xCr,2)
              U     = [U; element_disp(e,pos(:,k),enrich_node(:,k),u,k)];
          end
          % compute derivatives of u w.r.t xy
          H(1,1) = B(1,1:2:leB)*U(1:2:leB);    % u,x
          H(1,2) = B(2,2:2:leB)*U(1:2:leB);    % u,y
          H(2,1) = B(1,1:2:leB)*U(2:2:leB);    % v,x
          H(2,2) = B(2,2:2:leB)*U(2:2:leB);    % v,y
          
          
          epsilon = B*U ;
          sigma   = C*epsilon;
          
          voit2ind    = [1 3;3 2];
          gradqloc    = QT*gradq';
          graddisploc = QT*H*QT';
          stressloc   = QT*sigma(voit2ind)*QT';
          
          % ++++++++++++++++++
          %  Auxiliary fields
          % ++++++++++++++++++
          
          xp    = QT *(Gpt - xyTip)';           % local coordinates
          r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2)) ;
          theta = atan2(xp(2),xp(1)) ;
          
          for mode = 1:2

              [AuxStress,AuxGradDisp,AuxEps] = f_auxiliary(xp,r,theta,mu,kappa,mode);
              
              % +++++++++++++++
              %   J integral
              % +++++++++++++++
              I1= (stressloc(1,1) * AuxGradDisp(1,1) + stressloc(2,1) * AuxGradDisp(2,1) ) * gradqloc(1) + ...
                  (stressloc(1,2) * AuxGradDisp(1,1) + stressloc(2,2) * AuxGradDisp(2,1) ) * gradqloc(2);
              
              I2= (AuxStress(1,1) * graddisploc(1,1) + AuxStress(2,1) * graddisploc(2,1) ) * gradqloc(1) + ...
                  (AuxStress(2,1) * graddisploc(1,1) + AuxStress(2,2) * graddisploc(2,1) ) * gradqloc(2);
              
              StrainEnergy = 0;
              for i=1:2 %size(AuxEpsm1,1)
                  for j=1:2  %size(AuxEpsm1,2)
                      StrainEnergy = StrainEnergy +  stressloc(i,j)*AuxEps(i,j);
                  end
              end
              
              % Interaction integral I
              I(mode,1) = I(mode,1) + (I1 + I2 - StrainEnergy*gradqloc(1))*det(J0)*wt;
          end   %loop on mode
      end       % of quadrature loop
    end

    if wall_force
      if melange
        if ismember(e,mel_elems)
          mel_bool = 1;
        end
      end
      if ( ismember(e,vertex_elem) || ismember(e,split_elem) || ismember(e,tip_elem) )  && (any(elem_force(1,e,:)) | mel_bool)
        % The I integral needs to be adjusted to account for forces on the rift wall
        [ap,apg] = f_crack_wall(e,nnode,corner,tip_elem,vertex_elem,elem_crk,xyTip,xVertex,crack_nodes); % elem_crk in natural coordinates
        ap = f_align_lp_gc(ap,[apg(1,:),apg(end,:)],sctr);
        %keyboard
        for seg = 1:length(ap) - 1
          % find the distance between the two intersects (should be able to do this with det(J)
          p = ap(seg:seg+1,:);
          pg = [apg(seg,:),apg(seg+1,:)];
          [W,Q] = quadrature(wall_int,'GAUSS',1) ;
          % find the distance between the two intersects (should be able to do this with det(J)
          [l,nv,mv,nnt,nmt,mmt] = f_segment_dist(pg);
          %if flag_end == 1
            %nv = -1*nv;
          %end

          JO = l/2;
          nlocal = QT*nv';

          if mel_bool
            ii = find(mel_elems(:,2)==e);
            mT = mel_elems(ii,3);
            kn = mel_elems(ii,1); 
            JN = [ mv',nv' ]; % this is the rotation matrix
            [A,~,~,~,~] = f_enrich_assembly(e,pos,type_elem,elem_crk,enrich_node);
            U = u(A);
            B = [ ] ;
            Gpt = [1/3,1/3]; %CST 
            [N,dNdxi] = lagrange_basis(elemType,Gpt) ;
            [B, dJ] = xfemBmel(Gpt,e,type_elem,enrich_node(:,1),elem_crk,xVertex,xyTip,crack_nodes,kn,JN,mT,mE);
            eps_sub = B*U ;
            mstress = Cm1*eps_sub ;
            mel_stress = [mstress(1) mstress(3); mstress(3) mstress(2)];
            mel_local = QT*mel_stress*QT';
          end



          % ----------------------------
          % start loop over Gauss points
          % -----------------------------
          for gq = 1:size(W,1)
            [N1,dNdx1]=lagrange_basis('L2',Q(gq));
            pt = N1'*p;
            wt = W(gq,:);
            [N,dNdxi] = lagrange_basis(elemType,pt);
            Gpt = N' * node(sctr,:);     % GP in global coord
            % ++++++++++++++
            % q at pt 
            % ++++++++++++++
            qm2 = N'*qf';
            %qm1 = nv*qpt;
            %qm2 = -1*nv*qpt;

            % ++++++++++++++++++
            % stress at crack lip 
            % ++++++++++++++++++
            sig_elem = [0 elem_force(seg,e,2*gq); elem_force(seg,e,2*gq) elem_force(seg,e,2*gq-1)];
            %angl = atan2(nv(2),nv(1));
            %rotQT = [ cos(angl), sin(angl); -sin(angl), cos(angl) ]; 
            %sig_global = rotQT*sig_elem*rotQT';
            %sig_local = QT*sig_global*QT';
            sig_local1 = sig_elem/-2;
            sig_local2 = sig_elem/2;
            if mel_bool
              sig_local1 = sig_elem/-2 - mel_local;
              sig_local2 = sig_elem/2 + mel_local;
            end

            % ++++++++++++++++++
            %  Auxiliary fields
            % ++++++++++++++++++
            
            xp    = QT *(Gpt - xyTip)';           % local coordinate to tip
            r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2)) ;
            theta = pi;
            % if theta is equal to pi then we have to force one side to be consistently positive
            %if abs(abs(theta)-pi) < 1e-8 
              %theta = pi;
            %end

            for mode = 1:2
                [AuxStress,AuxGradDisp,AuxEps] = f_auxiliary(xp,r,theta,mu,kappa,mode);
            
                % +++++++++++++++
                %  Surface part of the I integral 
                % +++++++++++++++
                %keyboard
                I_wall1 = (sig_local1(1,2) * AuxGradDisp(1,1) + sig_local1(2,2) * AuxGradDisp(2,1) ) * qm2;
                %keyboard
                % Interaction integral I
                %keyboard
                theta =-1*pi;
                [AuxStress,AuxGradDisp,AuxEps] = f_auxiliary(xp,r,theta,mu,kappa,mode);

                % +++++++++++++++
                %  Surface part of the I integral 
                % +++++++++++++++
                I_wall2= (sig_local2(1,2) * AuxGradDisp(1,1) + sig_local2(2,2) * AuxGradDisp(2,1) ) * qm2;
                %keyboard
                
                % Interaction integral I
                Iw(mode,1) = Iw(mode,1) + I_wall1*det(JO)*wt + I_wall2*det(JO)*wt;
            end   %loop on mode
          end       % of quadrature loop
        end %segments in element
      end % if split_node and there is a force on wall 
    end
    
    for ni = 1:nn
      nod = sctr(ni);
      Gpt = node(nod,:);
      % ++++++++++++++++++
      % q at nodes 
      % ++++++++++++++++++
      pt = nnode(ni,:);  
      [N,dNdxi] = lagrange_basis(elemType,pt);
      qn = N'*qf';

      if ~ismember(nod,done_nodes) 
        done_nodes = [nod, done_nodes];



        % ++++++++++++++++++
        %  Auxiliary fields
        % ++++++++++++++++++
        xp    = QT *(Gpt - xyTip)';           % global coordinate to tip
        r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2)) ;
        theta = atan2(xp(2),xp(1)) ;

        for mode = 1:2
          [AuxStress,AuxGradDisp,AuxEps] = f_auxiliary(xp,r,theta,mu,kappa,mode);
          %pe = plot(Gpt(1),Gpt(2),'y*','markersize',10);
          f_i = QT * F(2*nod-1:2*nod);
          Fdudx = (f_i(1) * AuxGradDisp(1,1) + f_i(2)*AuxGradDisp(2,1))*qn;
          If(mode,1) = If(mode,1) - Fdudx;
        end
      end
    end
end
Knum = I.*E / (2*(1-nu^2)) % plain strain
KI = Knum(1);
KII = Knum(2);
theta_inc = 2 * atand((-2*KII / KI) / (1 + sqrt(1 + 8 * (KII / KI) ^ 2))); %deg
theta_inc = theta_inc * pi / 180.;%rad

kstr = ['Tip',num2str(flag_end),': stress-strain contribution K1 is ',num2str(KI),'   K2 is ',num2str(KII),'  and theta is ',num2str(theta_inc),'\n'];
fprintf(output_file,kstr)

Knumf = If.*E / (2*(1-nu^2)) % plain strain
KI = Knumf(1);
KII = Knumf(2);
theta_inc = 2 * atand((-2*KII / KI) / (1 + sqrt(1 + 8 * (KII / KI) ^ 2))); %deg
theta_inc = theta_inc * pi / 180.;%rad
kstr = ['Tip',num2str(flag_end),':  body forces contribution K1 is ',num2str(KI),'   K2 is ',num2str(KII),'  and theta is ',num2str(theta_inc),'\n'];
fprintf(output_file,kstr)

if any(Iw)
  Knumw = Iw.*E / (2*(1-nu^2)) % plain strain
  KI = Knumw(1);
  KII = Knumw(2);
  theta_inc = 2 * atand((-2*KII / KI) / (1 + sqrt(1 + 8 * (KII / KI) ^ 2))); %deg
  theta_inc = theta_inc * pi / 180.;%rad
  kstr = ['Tip',num2str(flag_end),': wall forces contribution K1 is ',num2str(KI),'   K2 is ',num2str(KII),'  and theta is ',num2str(theta_inc),'\n'];
  fprintf(output_file,kstr)
else
  Knumw = zeros(2,1);
end

Knum = Knum + Knumf + Knumw;
KI = Knum(1);
KII = Knum(2);
theta_inc = 2 * atand((-2*KII / KI) / (1 + sqrt(1 + 8 * (KII / KI) ^ 2))); %deg
theta_inc = theta_inc * pi / 180.;%rad
kstr = ['Tip',num2str(flag_end),': K1 is ',num2str(KI),'   K2 is ',num2str(KII),'  and theta is ',num2str(theta_inc),'\n'];
fprintf(output_file,kstr)



