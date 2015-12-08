%
% just test whether x^3+x^2 can be optimized over 
%

function stat = ipopt_test

%
% You need to update the path to IPOPT!
%
warning('off','MATLAB:dispatcher:pathWarning')
switch computer 
 case 'GLNX86' 
  addpath('./ipopt/ipopt-3.3.5-mumps-i686/matlab');
case 'GLNXA64' % cluster
  addpath('./ipopt/ipopt-3.3.5-mumps-amd64/matlab');
end  
warning('on','MATLAB:dispatcher:pathWarning')

x0 = 7;
lb = 2; 
ub = 10;
lbc = [];
ubc = [];

tolerance = 1e-4;

% call MKL with fininte number of kernels Ks
[x,multiplier,numiter] = ipopt(x0,lb,ub,lbc,ubc,...
			       @(x) compute_objective(x),...
			       @(x) compute_gradient(x),...
			       @(x) compute_constraints(x),...
			       @(x, rso) compute_jacobian(x, rso),...
			       '', [],'',[],...
			       'jac_c_constant','yes',...
			       'hessian_approximation','limited-memory',...
			       'mu_strategy','adaptive',...
			       'tol',tolerance,...
			       'print_level',0,... % 5 for verbosity
			       'max_iter',1000);

fprintf('solution = %g\n',x);
assert( abs((x-2))<10*tolerance);

fprintf('Test Passed\n');

return;


function f = compute_objective(x);
f = x.^3 + x.^2;

% gradient function
function df = compute_gradient(x)
df = 3*x.^2 + 2*x;

% constraints sum(w) = 1
function f = compute_constraints(x)
f = x;

% Jacobian
function J = compute_jacobian(x,return_structure_only);
J = sparse(ones(size(x)))';


