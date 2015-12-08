% Subproblem solver for Gaussian kernels with non-negative
% isotropic covariance matrix
%
% This functions aims to find local maxima of the function
%	T(\theta) = \sum_i\sum_j \alpha_i\alpha_j k(x_i,x_j;\theta)
% where \alpha can be positive or negative. For this function the
% kernel function is of the form
%
% k(x,x';\theta) = \exp(- \sum_{d=1}^D \theta_d (x_d-x'_d)^2 )
%
%
% Usage: 
%	[K,new_facs,parm] = IKLsubproblem_nonnegisotropic(operation,alpha,lambda,Ds,parm,w,parm_range,opts)
%
% Input: (needed for IKLmain!)
%	operation - either 'init' : generate the starting matrices
%			'search' : search for violating constraint
%			'build' : build a Gram matrix given the data
%	alpha - current SVM multipliers alpha
%	lambda - Lagrange multiplier of the constraint sum(w)=1,
%	  this is the threshold for violating constraints (T(theta)>lambda)
%	w - mixing coefficients (only using 'build')
%	parm - subproblem parameters, startpoints, points included,
%		etc. (can differ for any subproblem)
%
% Input: (user defined, specific to the subproblem)
%	Ds - distances between training data
%	parm_range - admissible range for kernel parameters
%	opts - opts.eps (1e-3) = Tolerance for Ipopt solver (for subproblem)
%
% Output: 
%	K - generated kernel matrices [n x n x M] or [n x n_test]
%		for 'build'
%	new_facs - new kernel parameters
%	parm - additional data untouched from IKL
%
%
% Peter Gehler 07/2008 pgehler@tuebingen.mpg.de

function [K,new_facs,parm] = IKLsubproblem_nonnegisotropic(operation,al,lam,Ds,parm,w,parm_range,opts)

if ~exist('opts','var'), opts = []; end

%Multiple Pricing: how many constraints to return in one search
if ~isfield(opts,'constraints_per_iteration')
    opts.constraints_per_iteration = 1;
end

% Tolerance for Ipopt solver 
if ~isfield(opts,'eps')
    opts.eps = 1e-3;
end

% The violation must be exceed this threshold in order to include
% the corresponding constraint.
if ~isfield(opts,'constraint_eps')
    opts.constraint_eps = 1e-2; 
end

warning('off','MATLAB:dispatcher:pathWarning')
addpath('~/software/libsvm-matlab');
switch computer 
 case 'GLNX86' 
  addpath('./ipopt/ipopt-3.3.5-mumps-i686/matlab');
case 'GLNXA64' % cluster
  addpath('./ipopt/ipopt-3.3.5-mumps-amd64/matlab');
end  
warning('on','MATLAB:dispatcher:pathWarning')


DO_TEST = 1;

if ~exist('parm_range','var') & ~(strcmp(operation,'build')>0)
    warning('no range for valid exponents given: setting to [0.1,4]');
    parm_range = [0.1;4];
end

nFactors = size(Ds,3);
new_facs = [];


switch operation
 case 'init'
  
  %
  % initialization : 
  % all single kernel matrices
  % + product of them
  %
  % 3Dimensional example with parm_range [1,4], the start points will be
  %
  % [1,0,0]
  % [2,0,0]
  % [3,0,0]
  % [4,0,0]
  % [0,1,0]
  % ...
  % [0,0,4]
  % [1,1,1]
  % [2,2,2]
  % [3,3,3]
  % [4,4,4]

  cntr = 1;
  %parm.facs = zeros(size(Ds,3),size(Ds,3));

  sig0 = mean(parm_range);

  if nFactors ==1
      scales = linspace(parm_range(1),parm_range(2),50);
  else
      scales = linspace(parm_range(1),parm_range(2),4);
  end	
  start_points = zeros(numel(scales)*nFactors,nFactors);
  
  K = [];
  for k=1:nFactors
      %tD = Ds(:,:,k);
      %K(:,:,k) = exp(-sig0 * Ds(:,:,k));
      %parm.facs(k,k) = sig0;
      for kk=1:numel(scales)
	  start_points(cntr,k) = scales(kk);
	  cntr = cntr + 1;
      end
  end

  
  tD = sum(Ds,3)/nFactors;
  sig0 = mean(parm_range);
  %K(:,:,end+1) = exp(-sig0*tD);
  %parm.facs(end+1,1:nFactors) = sig0;

  K(:,:,1) = exp(-sig0*tD);
  parm.facs(1,1:nFactors) = sig0;

  
  if nFactors > 1
      for kk=1:numel(scales)
	  start_points(cntr,:) = scales(kk);
	  cntr = cntr + 1;
      end
  end

  start_points = removeDuplicates(start_points')';
  
  new_facs = parm.facs;
  parm.startpoints = start_points;

  fprintf('SP: generated %d Gram matrices\n',size(K,3));
  fprintf('SP: generated %d start_points\n',size(start_points,1));
  
  assert(size(start_points,2)==nFactors);

 case 'build'

  % PREDICTION
  % build the Gram matrix and return
  assert(ndims(Ds)==2 | (size(Ds,3)==size(parm,2)),'dim mismatch');
  K = 0;
  for k=1:size(parm,1)
      d = 0;
      for j=1:size(parm,2)
	  d = d + parm(k,j) * Ds(:,:,j);
      end
      K = K + w(k) * exp(-d);
  end

 case 'search'
  

  K = []; % matrix to return

  svset = find(al~=0);
  al = al(svset);


  % lower and upper bounds for exponents
  lb = parm_range(1) * ones(nFactors,1);
  ub = parm_range(2) * ones(nFactors,1);


  % specify the starting points ...
  if isfield(parm,'startpoints')
      startpoints = [parm.startpoints;parm.facs]';
  else
      startpoints = parm.facs'; 
  end
  startpoints = removeDuplicates(startpoints);
  nof_startpoints = size(startpoints,2);
  shuff = randperm(nof_startpoints);
  startpoints = startpoints(:,shuff);

  % some global variabels for ipopt optimization
  global subprob_obj_counter subprob_deriv_counter subprob_hessian_counter
  global subprob_recompute_counter
  global subprob_f subprob_df subprob_fac subprob_H subprob_sig


  nof_violating_constraints = 0;
  for n=1:nof_startpoints
      fac0 = startpoints(:,n);
      
      subprob_obj_counter=0;
      subprob_deriv_counter=0;
      subprob_hessian_counter=0;
      subprob_recompute_counter = 0;
      
      subprob_fac = Inf;
      subprob_f = -Inf;
      subprob_df = -Inf;
      subprob_H = -Inf;
      subprob_sig = Inf;

      [fac0,multiplier,numiter] = ipopt(fac0,lb,ub,[],[],...
					@(fac) compute_objective(fac,al,Ds,svset),...
					@(fac) compute_gradient(fac,al,Ds,svset),...
					@(fac) compute_constraints(fac,al,Ds,svset),...
					@(fac, rso) compute_jacobian(fac, rso),...
					@(fac, sigma, lambda,rso) compute_hessian(fac, sigma, lambda, rso, al,Ds,svset),...
					[],'',[],...
					'hessian_approximation','exact',...
					'jac_c_constant','yes',...
					'mu_strategy','adaptive',...
					'tol',opts.eps,...
					'print_level',0,...
					'max_iter',20);
      numiter = numiter + 1; %C->Matlab
      assert(numiter == subprob_deriv_counter,'#derivative calls != #ipopt  iterations');

      
      fprintf('SP: IpOpt iter %d, ',numiter);
      fprintf('Evaluations: F:%d, dF:%d, H:%d, ',subprob_obj_counter,subprob_deriv_counter,subprob_hessian_counter);
      fprintf('Computations needed: %d\n',subprob_recompute_counter);

      if (DO_TEST)
	  tmp_f = subprob_f;
	  subprob_fac = -Inf(nFactors,1);
	  f = compute_objective(fac0,al,Ds,svset);
	  assert(tmp_f==f|abs(tmp_f-f)./f<1e-5,'objective differs');
      end

      % check if the found local maxima violates the contraint
      % pg : to do: criterion relative to absolute value
      if (-subprob_f - lam) >= opts.constraint_eps
	  
	  % yes it did ...
	  fprintf('found violating constraint...');

	  % ... already in the list? 
	  already_included = any(dist_euclid(parm.facs,fac0')<opts.eps);

	  if already_included
	      fprintf('already included: skipping\n');
	      continue;
	  end
	  
	  % ... add it to the constraint set
	  fprintf('adding it to the constraint set\n');
	  fprintf('next factors');    fprintf(' %.3f',fac0);fprintf('\n');
	  
	  D = 0;
	  for k=1:nFactors
	      D = D + fac0(k) * Ds(:,:,k);
	  end
	  K(:,:,nof_violating_constraints+1) = exp(-D);
	  
	  parm.facs(end+1,:) = fac0;
	  new_facs(end+1,:) = fac0;
	  
	  nof_violating_constraints = nof_violating_constraints + 1;
	  
	  % ... are there enough constraints found ? ...
	  if (nof_violating_constraints >= opts.constraints_per_iteration)
	      fprintf('SP: found maximum of %d constraints, returning\n',nof_violating_constraints);
	      return;
	  end
      end
  end


  % ... no violating constraint has been found, return nothing ...
  if nof_violating_constraints == 0
      K = [];
      new_facs = [];
  else
      fprintf('SP: found only %d constraints (of %d requested) returning\n',nof_violating_constraints,opts.constraints_per_iteration);
  end

  return;
end



% Objective
function f = compute_objective(fac,al,Ks,svset)
global subprob_obj_counter
subprob_obj_counter = subprob_obj_counter + 1;
f = update_objective(fac,al,Ks,svset);


%Gradient
function df = compute_gradient(fac,al,Ks,svset)
global subprob_deriv_counter 
subprob_deriv_counter = subprob_deriv_counter +1;

[foo,df] = update_objective(fac,al,Ks,svset);



% Hessian
function H = compute_hessian(fac,ipopt_sigma,ipopt_lambda,return_structure_only,al,D,svset)
global subprob_hessian_counter
subprob_hessian_counter = subprob_hessian_counter +1;

if return_structure_only
    H = sparse(tril(ones(numel(fac),numel(fac))));
else
    [foo,foo2,H] = update_objective(fac,al,D,svset);
    H = ipopt_sigma * sparse(H);
end


% Constraints
function g = compute_constraints(fac)
g=fac;

% Jacobian
function J = compute_jacobian(fac,return_structure_only);
J = sparse(numel(fac),1);

%
% compute the objective, the gradient 
%
function [f,df,H] = update_objective(fac,al,Ds,svset)

assert(numel(fac)==size(Ds,3),'number of arguments != number of distances');
%assert(all(fac>=-1e-5),'factor turned negative');
%fac = max(fac,0);

global subprob_f subprob_df subprob_fac subprob_H

%
% All what is following should be computed in a mex file for
% more speed
% 

if any(size(fac)~=size(subprob_fac)) | any(fac~=subprob_fac)

    global subprob_recompute_counter
    subprob_recompute_counter = subprob_recompute_counter + 1;
    
    [subprob_f,subprob_df,subprob_H] = subproblem_product_objective(double(Ds(svset,svset,:)),fac,al);
    
    subprob_fac = fac;
end


f = subprob_f;
if nargout>1
    df = subprob_df;
    
    if nargout >2
	H = subprob_H;
    end
end


% remove the duplicate rows of the matrix points e.g. 
% 
% [1 0 0   => [1 0
%  0 0 0]      0 0]
function points = removeDuplicates(points);

D = dist_euclid(points');
[i1,i2] = ind2sub(size(D),find(D<1e-6));

indx = find(i1<=i2);
i1(indx) = [];
i2(indx) = [];

indx = [];
for i=1:max(i2)
    indx = [indx;i1(find(i2==i))];
end
indx = unique(indx);

points(:,indx) = [];

% sanity check
D = dist_euclid(points');
assert(sum(D(:)<1e-6)==size(D,1));

