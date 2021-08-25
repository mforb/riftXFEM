function roots = quadraticSolver(a,b,c)
  % quadratic solver to understand tests

  if ~isa(a,'numeric') || ~isa(b,'numeric') || ~isa(c,'numeric')
    error('quadraticSolver:InputMustBeNumeric', ...
         'Coefficients must be numeric.');
  end

  roots(1) = ( -b + sqrt( b*b - 4*a*c) )/ 2*a;
  roots(2) = ( -b - sqrt( b*b - 4*a*c) )/ 2*a;
end

