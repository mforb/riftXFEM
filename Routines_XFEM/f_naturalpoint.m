function [ nt ] = f_naturalpoint( pt, nodes, iter, epsilon )
% This MATLAB function was created by Martin Forbes (martin.forbes@postgrad.otago.ac.nz)
% The date of creation: Fri Nov 26 11:33:23 NZDT 2021
global elemType

inc = 1;
if strcmp(elemType,'Q4')
  coord = zeros(1,2);
elseif strcmp(elemType,'T3')
  coord = ones(1,2)*1/3;
else
  error('Element type not implemented (finding natural coord vertex/tip)') 
end
ksi   = 0;
eta   = 0;
converge = 0;

while (inc < iter)
    [N,dNdxi]=lagrange_basis(elemType,coord);   % compute shape functions
    x = N'*nodes(:,1);
    y = N'*nodes(:,2);
    df1dr = dNdxi(:,1)' * nodes(:,1);
    df1ds = dNdxi(:,2)' * nodes(:,1);
    df2dr = dNdxi(:,1)' * nodes(:,2);
    df2ds = dNdxi(:,2)' * nodes(:,2);
 
    f1 = x - pt(1);
    f2 = y - pt(2);

    detF = df1dr*df2ds - df1ds*df2dr ;  %Jacobien

    invf(1,1) =  1.0/detF * df2ds;
    invf(1,2) = -1.0/detF * df1ds;
    invf(2,1) = -1.0/detF * df2dr;
    invf(2,2) =  1.0/detF * df1dr;

    ksi = ksi - invf(1,1)*f1 - invf(1,2)*f2;
    eta = eta - invf(2,1)*f1 - invf(2,2)*f2;


    if( (abs(ksi - coord(1)) < epsilon) && ...
            (abs(eta - coord(2)) < epsilon) )
        inc  = iter + 1;
        nt = coord;
        converge = 1;
    else
        coord(1) = ksi;
        coord(2) = eta;
        inc = inc + 1;
    end
end
if ~converge
  warning('Natural point location did not converge')
end
end
