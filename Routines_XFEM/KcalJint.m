function [Knumerical,ThetaInc,xCr,stop] = KcalJint(xCr,...
    type_elem,enrdomain,elem_crk,enrich_node,crack_nodes,xVertex,...
    vertex_elem,pos,u,F,ipas,delta_inc,Knumerical,ThetaInc,...
    tip_elem,split_elem,corner_elem, tan_elem,elem_force,gn_inters)

global node element elemType
global E nu C sigmato
global Jint iMethod
global output_file
global quick_freeze
global typeProblem

stp1 = 0;
stp2 = 0;

% ThetaInc = [] ;
kstr = ['Calculating SIFs and propagating rift\n'];
for kk = 1:size(xCr,2) %what's the crack?
    disp([num2str(toc),'      Crack n. ',num2str(kk)])
    if ~isempty(Knumerical) ;
      K1_num = Knumerical{kk,1};
      K2_num = Knumerical{kk,2};
    else
      K1_num = [];
      K2_num = [];
    end
    if ~isempty(ThetaInc) ;
      ti1 = ThetaInc{kk,1};
      ti2 = ThetaInc{kk,2};
    else
      ti1 = [];
      ti2 = [];
    end

    xCr(kk).coornew1 = [];
    xCr(kk).coornew2 = [];
    tip = find(type_elem(:,kk) == 1);
    split = find(type_elem(:,kk) == 2);
    %find out the element containing the...
    for ii=1:size(tip_elem,1)
        iel = tip_elem(ii) ;
        sctr=element(iel,:);
        vv = node(sctr,:);
        flag1 = inhull(xCr(kk).coor(1,:),vv,[],1e-8); % left tip!
        if flag1 == 1
            seg   = xCr(kk).coor(1,:) - xCr(kk).coor(2,:);
            alpha = atan2(seg(2),seg(1));
            if strcmp(typeProblem,'ISSM') || strcmp(typeProblem,'centreF')
              [Knum,theta_inc] = SIF_BF(C,1,iel,elem_crk,xCr,type_elem,...
                  enrich_node,crack_nodes,xVertex,pos,u,F,kk,alpha,tip_elem,split_elem,vertex_elem,corner_elem,tan_elem,elem_force) ;
            else
              [Knum,theta_inc] = SIF(C,1,iel,elem_crk,xCr,type_elem,...
                  enrich_node,crack_nodes,xVertex,pos,u,kk,alpha,tip_elem,split_elem,vertex_elem,corner_elem,tan_elem,elem_force) ;
            end

            theta_inc = -1*theta_inc; 
            K1_num = [K1_num, Knum] ;
            ti1 = [ti1, theta_inc] ; % because of flipped y-axis 
            kstr = ['Tip 1 modified: K1 is ',num2str(Knum(1)),'   K2 is ',num2str(Knum(2)),'  and theta is ',num2str(theta_inc),'\n'];
            if xCr(kk).tip(1) & Knum(1) > 0
              [xCr(kk),theta_inc,kstr] = f_advance_tip(xCr(kk),1,delta_inc,theta_inc,alpha,kstr);
            else
              stp1 = 1;
            end
            fprintf(output_file,kstr)
        end
        flag2 = inhull(xCr(kk).coor(size(xCr(kk).coor,1),:),vv,[],1e-8); % and the right tip!
        if flag2 == 1
            seg   = xCr(kk).coor(size(xCr(kk).coor,1),:) - xCr(kk).coor(size(xCr(kk).coor,1)-1,:);
            alpha = atan2(seg(2),seg(1)) ;

            if strcmp(typeProblem,'ISSM') || strcmp(typeProblem,'centreF')
              [Knum,theta_inc] = SIF_BF(C,2,iel,elem_crk,xCr,type_elem,...
                  enrich_node,crack_nodes,xVertex,pos,u,F,kk,alpha,tip_elem,split_elem,vertex_elem,corner_elem,tan_elem,elem_force) ;
            else
              [Knum,theta_inc] = SIF(C,2,iel,elem_crk,xCr,type_elem,...
                  enrich_node,crack_nodes,xVertex,pos,u,kk,alpha,tip_elem,split_elem,vertex_elem,corner_elem,tan_elem,elem_force) ;
            end
            K2_num = [K2_num, Knum] ;
            if xCr(kk).tip(2) & Knum(1) > 0
              [xCr(kk),theta_inc,kstr] = f_advance_tip(xCr(kk),2,delta_inc,theta_inc,alpha,kstr);
            else
              stp2 = 1;
            end
            ti2 = [ti2, theta_inc] ;
            fprintf(output_file,kstr)
        end
    end
    xCr(kk).coor = [xCr(kk).coornew1;xCr(kk).coor;xCr(kk).coornew2] ;
    if isfield(xCr,'melange')
      mn1 = []; mn2 = [];
      w1 = []; w2 = [];
      if ~isempty(xCr(kk).coornew1)
        if quick_freeze
          mn1 = 1 ;
          w1 = 0;
        else
          mn1 = 0 ;
          w1 = 0;
        end
      end
      if ~isempty(xCr(kk).coornew2)
        if quick_freeze
          mn2 = 1;
          w2 = 0 ;
        else
          mn2 = 0;
          w2 =  0 ;
        end
      end
      xCr(kk).width = xCr(kk).width + gn_inters;
      xCr(kk).melange = [ mn1; xCr(kk).melange; mn2 ];
      xCr(kk).width = [ w1 , xCr(kk).width, w2 ];
      %if quick_freeze
        %xCr(kk).width = xCr(kk).width + 10*ones(size(xCr(kk).width)) ;
      %end 
    end
      

    ThetaInc{kk,1} = ti1;
    ThetaInc{kk,2} = ti2;
    Knumerical{kk,1} = K1_num;
    Knumerical{kk,2} = K2_num;
    stop = stp1*stp2;
end %kk
