function [Knum,theta_inc] = SIF(C,flag_end,tip,elem_crk,xCr,type_elem,enrich_node,crack_nodes,xVertex,pos,u,kk,alpha,...
    tip_elem,split_elem,vertex_elem,corner_elem,elem_force)

global node element elemType E nu
global iMethod iParam
global incR xc yc phiN
global lambda1 lambda2 nu1 nu2
global elemType typeMesh typeProblem typeCrack stressState
global wall_int

plotNode = 'NO' ;

if strcmp(elemType,'Q4') 
  intType = 'GAUSS' ;
  corner = [1 2 3 4 1] ;
  nnode = [-1 -1;1 -1;1 1;-1 1] ;
else
  intType = 'TRIANGULAR';
  corner = [1 2 3 1] ;
  nnode = [0 0;1 0;0 1] ;
end

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
elseif flag_end == 2
  xyTip = [elem_crk(tip,3) elem_crk(tip,4)] ;
end

[Jdomain,JWdomain,qnode,qnode2,radius] = Jdomainf(tip,xyTip,enrich_node);

I1 = 0;
I2 = 0;
I  = [zeros(2,1)];

QT  =[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];           % for the transformation to local coordinate

% ---------------------------
% Starting LOOP over ELEMENTS
%----------------------------
compt=0;
q=[];

for iel = 1 : size(Jdomain,2)
    e = Jdomain(iel) ; % current element
    sctr = element(e,:);
    nn   = length(sctr);
    normal_order = 6;          % max = 8
    tip_order    = 20;         % is not used -> if >20 allow a error flag (bug) if the J integrale is performed on blending's
    split_order  = 7;
    vertex_order = 7;
    
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
        
        
        % ++++++++++++++
        % Gradient of q
        % ++++++++++++++
        q     = qnode(iel,:);
        gradq = q*dNdx;
        
        % ++++++++++++++
        % Stress at GPs
        % ++++++++++++++
        
%                 % Compliance matrix
%         x = pt(1);
%         y = pt(2);
%         d = sqrt((x-xc)^2+(y-yc)^2);
%         if d < incR 
%             E  = lambda2 ;
%             nu = nu2;
%         else
%             E  = lambda1 ;
%             nu = nu1;
%         end
%         
%         if ( strcmp(stressState,'PLANE_STRESS') )
%             C=E/(1-nu^2)*[ 1   nu 0;
%                 nu  1  0 ;
%                 0   0  0.5*(1-nu) ];
%         else
%             C=E/(1+nu)/(1-2*nu)*[ 1-nu  nu  0;
%                 nu    1-nu 0;
%                 0     0  0.5-nu ];
%         end
        
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

        K1 = 1.0 ;
        K2 = K1  ;
        
        mu = E/(2.+ nu + nu);
        kappa = 3-4*nu;    %Kolosov coeff, Plain strain
        
        SQR  = sqrt(r);
        CT   = cos(theta);
        ST   = sin(theta);
        CT2  = cos(theta/2);
        ST2  = sin(theta/2);
        C3T2 = cos(3*theta/2);
        S3T2 = sin(3*theta/2);
        
        drdx = CT;
        drdy = ST;
        dtdx = -ST/r;
        dtdy = CT/r;
        
        FACStress1 = sqrt(1/(2*pi));
        FACStress2 = FACStress1;
        
        FACDisp1 = sqrt(1/(2*pi))/(2*mu);
        FACDisp2 = FACDisp1;
        
        AuxStress   = zeros(2,2);
        AuxGradDisp = zeros(2,2);
        AuxEps      = zeros(2,2);
        
        for mode = 1:2
            if (mode == 1)
                
                AuxStress(1,1) = K1*FACStress1/SQR*CT2*(1-ST2*S3T2);
                AuxStress(2,2) = K1*FACStress1/SQR*CT2*(1+ST2*S3T2);
                AuxStress(1,2) = K1*FACStress1/SQR*ST2*CT2*C3T2;
                AuxStress(2,1) = AuxStress(1,2);
                
                u1    = K1*FACDisp1*SQR*CT2*(kappa - CT);
                du1dr = K1*FACDisp1*0.5/SQR*CT2*(kappa - CT);
                du1dt = K1*FACDisp1*SQR*(-0.5*ST2*(kappa - CT) + CT2*ST);
                
                u2    = K1*FACDisp1*SQR*ST2*(kappa - CT);
                du2dr = K1*FACDisp1*0.5/SQR*ST2*(kappa - CT);
                du2dt = K1*FACDisp1*SQR*(0.5*CT2*(kappa - CT) + ST2*ST);
                
                AuxGradDisp(1,1) = du1dr*drdx + du1dt*dtdx;
                AuxGradDisp(1,2) = du1dr*drdy + du1dt*dtdy;
                AuxGradDisp(2,1) = du2dr*drdx + du2dt*dtdx;
                AuxGradDisp(2,2) = du2dr*drdy + du2dt*dtdy;
                
                AuxEps(1,1) = AuxGradDisp(1,1);
                AuxEps(2,1) = 0.5*(AuxGradDisp(2,1) + AuxGradDisp(1,2));
                AuxEps(1,2) = AuxEps(2,1);
                AuxEps(2,2) = AuxGradDisp(2,2);
                
            elseif (mode == 2)
                AuxStress(1,1) = -K2*FACStress2/SQR*ST2*(2+CT2*C3T2);
                AuxStress(2,2) = K2*FACStress2/SQR*ST2*CT2*C3T2;
                AuxStress(1,2) = K2*FACStress2/SQR*CT2*(1-ST2*S3T2);
                AuxStress(2,1) = AuxStress(1,2);
                
                u1    = K2*FACDisp2*SQR*ST2*(kappa + 2 + CT);
                du1dr = K2*FACDisp2*0.5/SQR*ST2*(kappa + 2 + CT);
                du1dt = K2*FACDisp2*SQR*(0.5*CT2*(kappa + 2 + CT) - ST2*ST);
                
                u2    = -K2*FACDisp2*SQR*CT2*(kappa - 2 + CT);
                du2dr = -K2*FACDisp2*0.5*(1/SQR)*CT2*(kappa - 2 + CT);
                du2dt = -K2*FACDisp2*SQR*(-0.5*ST2*(kappa - 2 + CT) - CT2*ST);
                
                AuxGradDisp(1,1) = du1dr*drdx + du1dt*dtdx;
                AuxGradDisp(1,2) = du1dr*drdy + du1dt*dtdy;
                AuxGradDisp(2,1) = du2dr*drdx + du2dt*dtdx;
                AuxGradDisp(2,2) = du2dr*drdy + du2dt*dtdy;
                
                AuxEps(1,1) = AuxGradDisp(1,1);
                AuxEps(2,1) = 0.5*(AuxGradDisp(2,1) + AuxGradDisp(1,2));
                AuxEps(1,2) = AuxEps(2,1);
                AuxEps(2,2) = AuxGradDisp(2,2);
            end
            
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
Knum = I.*E / (2*(1-nu^2)) % plain strain

for iel = 1 : size(JWdomain,2)
    e = JWdomain(iel) ; % current element
    sctr = element(e,:);
    nn   = length(sctr);

    if ( ismember(e,vertex_elem) || ismember(e,split_elem) || ismember(e,tip_elem) )  && any(elem_force(1,e,:))
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

        JO = l/2;
        nlocal = QT*nv';



        % ++++++++++++++
        % q at nodes
        % ++++++++++++++
        q     = qnode2(iel,:);

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
          qpt = N'*q';
          qm1 = nv*qpt;
          qm2 = -1*nv*qpt;

          % ++++++++++++++++++
          % stress at crack lip 
          % ++++++++++++++++++
          sig_elem = [elem_force(seg,e,2*gq-1) elem_force(seg,e,2*gq); elem_force(seg,e,2*gq) 0];
          angl = atan2(nv(2),nv(1));
          rotQT = [ cos(angl), sin(angl); -sin(angl), cos(angl) ]; 
          sig_global = rotQT*sig_elem*rotQT';
          sig_local = QT*sig_global*QT';
          sig_local1 = sig_local/2;
          sig_local2 = sig_local/-2;

          % ++++++++++++++++++
          %  Auxiliary fields
          % ++++++++++++++++++
          
          xp    = QT *(Gpt - xyTip)';           % local coordinate to tip
          r     = sqrt(xp(1)*xp(1)+xp(2)*xp(2)) ;
          theta = atan2(xp(2),xp(1)) ;
          % if theta is equal to pi then we have to force one side to be consistently positive
          if abs(abs(theta)-pi) < 1e-8 
            theta = pi;
          end

          K1 = 1.0 ;
          K2 = K1  ;
          
          mu = E/(2.+ nu + nu);
          kappa = 3-4*nu;    %Kolosov coeff, Plain strain
          
          SQR  = sqrt(r);
          CT   = cos(theta);
          ST   = sin(theta);
          CT2  = cos(theta/2);
          ST2  = sin(theta/2);
          C3T2 = cos(3*theta/2);
          S3T2 = sin(3*theta/2);
          
          drdx = CT;
          drdy = ST;
          dtdx = -ST/r;
          dtdy = CT/r;
          for mode = 1:2
            if (mode == 1)
                  
                  u1    = K1*FACDisp1*SQR*CT2*(kappa - CT);
                  du1dr = K1*FACDisp1*0.5/SQR*CT2*(kappa - CT);
                  du1dt = K1*FACDisp1*SQR*(-0.5*ST2*(kappa - CT) + CT2*ST);
                  
                  u2    = K1*FACDisp1*SQR*ST2*(kappa - CT);
                  du2dr = K1*FACDisp1*0.5/SQR*ST2*(kappa - CT);
                  du2dt = K1*FACDisp1*SQR*(0.5*CT2*(kappa - CT) + ST2*ST);
                  
                  AuxGradDisp(1,1) = du1dr*drdx + du1dt*dtdx;
                  AuxGradDisp(1,2) = du1dr*drdy + du1dt*dtdy;
                  AuxGradDisp(2,1) = du2dr*drdx + du2dt*dtdx;
                  AuxGradDisp(2,2) = du2dr*drdy + du2dt*dtdy;
                  
                  AuxEps(1,1) = AuxGradDisp(1,1);
                  AuxEps(2,1) = 0.5*(AuxGradDisp(2,1) + AuxGradDisp(1,2));
                  AuxEps(1,2) = AuxEps(2,1);
                  AuxEps(2,2) = AuxGradDisp(2,2);
                  
              elseif (mode == 2)
                  
                  u1    = K2*FACDisp2*SQR*ST2*(kappa + 2 + CT);
                  du1dr = K2*FACDisp2*0.5/SQR*ST2*(kappa + 2 + CT);
                  du1dt = K2*FACDisp2*SQR*(0.5*CT2*(kappa + 2 + CT) - ST2*ST);
                  
                  u2    = -K2*FACDisp2*SQR*CT2*(kappa - 2 + CT);
                  du2dr = -K2*FACDisp2*0.5*(1/SQR)*CT2*(kappa - 2 + CT);
                  du2dt = -K2*FACDisp2*SQR*(-0.5*ST2*(kappa - 2 + CT) - CT2*ST);
                  
                  AuxGradDisp(1,1) = du1dr*drdx + du1dt*dtdx;
                  AuxGradDisp(1,2) = du1dr*drdy + du1dt*dtdy;
                  AuxGradDisp(2,1) = du2dr*drdx + du2dt*dtdx;
                  AuxGradDisp(2,2) = du2dr*drdy + du2dt*dtdy;
                  
                  AuxEps(1,1) = AuxGradDisp(1,1);
                  AuxEps(2,1) = 0.5*(AuxGradDisp(2,1) + AuxGradDisp(1,2));
                  AuxEps(1,2) = AuxEps(2,1);
                  AuxEps(2,2) = AuxGradDisp(2,2);
              end
              
              % +++++++++++++++
              %  Surface part of the I integral 
              % +++++++++++++++
              I_wall1 = -1*(sig_local1(1,2) * AuxGradDisp(1,1) + sig_local1(2,2) * AuxGradDisp(2,1) ) * qm1(2);
              % Interaction integral I
              %keyboard
              I(mode,1) = I(mode,1) + I_wall1*det(JO)*wt;
          end   %loop on mode

          theta =-1*(2*pi - theta);

          CT   = cos(theta);
          ST   = sin(theta);
          CT2  = cos(theta/2);
          ST2  = sin(theta/2);
          C3T2 = cos(3*theta/2);
          S3T2 = sin(3*theta/2);
          
          drdx = CT;
          drdy = ST;
          dtdx = -ST/r;
          dtdy = CT/r;
          for mode = 1:2
            if (mode == 1)
                  
                  u1    = K1*FACDisp1*SQR*CT2*(kappa - CT);
                  du1dr = K1*FACDisp1*0.5/SQR*CT2*(kappa - CT);
                  du1dt = K1*FACDisp1*SQR*(-0.5*ST2*(kappa - CT) + CT2*ST);
                  
                  u2    = K1*FACDisp1*SQR*ST2*(kappa - CT);
                  du2dr = K1*FACDisp1*0.5/SQR*ST2*(kappa - CT);
                  du2dt = K1*FACDisp1*SQR*(0.5*CT2*(kappa - CT) + ST2*ST);
                  
                  AuxGradDisp(1,1) = du1dr*drdx + du1dt*dtdx;
                  AuxGradDisp(1,2) = du1dr*drdy + du1dt*dtdy;
                  AuxGradDisp(2,1) = du2dr*drdx + du2dt*dtdx;
                  AuxGradDisp(2,2) = du2dr*drdy + du2dt*dtdy;
                  
                  AuxEps(1,1) = AuxGradDisp(1,1);
                  AuxEps(2,1) = 0.5*(AuxGradDisp(2,1) + AuxGradDisp(1,2));
                  AuxEps(1,2) = AuxEps(2,1);
                  AuxEps(2,2) = AuxGradDisp(2,2);
                  
              elseif (mode == 2)
                  
                  u1    = K2*FACDisp2*SQR*ST2*(kappa + 2 + CT);
                  du1dr = K2*FACDisp2*0.5/SQR*ST2*(kappa + 2 + CT);
                  du1dt = K2*FACDisp2*SQR*(0.5*CT2*(kappa + 2 + CT) - ST2*ST);
                  
                  u2    = -K2*FACDisp2*SQR*CT2*(kappa - 2 + CT);
                  du2dr = -K2*FACDisp2*0.5*(1/SQR)*CT2*(kappa - 2 + CT);
                  du2dt = -K2*FACDisp2*SQR*(-0.5*ST2*(kappa - 2 + CT) - CT2*ST);
                  
                  AuxGradDisp(1,1) = du1dr*drdx + du1dt*dtdx;
                  AuxGradDisp(1,2) = du1dr*drdy + du1dt*dtdy;
                  AuxGradDisp(2,1) = du2dr*drdx + du2dt*dtdx;
                  AuxGradDisp(2,2) = du2dr*drdy + du2dt*dtdy;
                  
                  AuxEps(1,1) = AuxGradDisp(1,1);
                  AuxEps(2,1) = 0.5*(AuxGradDisp(2,1) + AuxGradDisp(1,2));
                  AuxEps(1,2) = AuxEps(2,1);
                  AuxEps(2,2) = AuxGradDisp(2,2);
              end
              
              % +++++++++++++++
              %  Surface part of the I integral 
              % +++++++++++++++
              I_wall2= (sig_local2(1,2) * AuxGradDisp(1,1) + sig_local2(2,2) * AuxGradDisp(2,1) ) * qm2(2);
              %keyboard
              
              % Interaction integral I
              I(mode,1) = I(mode,1) + I_wall2*det(JO)*wt;
          end   %loop on mode
        end       % of quadrature loop
      end %segments in element
    end % if split_node and there is a force on wall 
end           % end of element loop



% Compute SIFs from I integral
Knum = I.*E / (2*(1-nu^2)) % plain strain
Knum = Roundoffa(Knum,5);
KI = Knum(1);
KII = Knum(2);
theta_inc = 2 * atand((-2*KII / KI) / (1 + sqrt(1 + 8 * (KII / KI) ^ 2))); %deg
theta_inc = theta_inc * pi / 180.;%rad


%
% figure
% hold on
% %plot the circle
% theta = -pi:0.1:pi;
% xo = xyTip(1) + radius*cos(theta) ;
% yo = xyTip(2) + radius*sin(theta) ;
% plot(xo,yo,'k-');
% %plot_mesh(node,element,elemType,'b-')
% %plot_mesh(node,element(Jdomain,:),elemType,'r-')
%
% plotMesh(node,element,elemType,'b-',plotNode);
% plotMesh(node,element(Jdomain,:),elemType,'b-',plotNode);
%
% for k = 1:size(xCr,2)
%   for kj = 1:size(xCr(k).coor,1)-1
%     cr = plot(xCr(k).coor(kj:kj+1,1),xCr(k).coor(kj:kj+1,2),'r-');
%     set(cr,'LineWidth',3);
%   end
% end
% pause
% % -------------------------------------
% % -------------------------------------
