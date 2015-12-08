% Infinite Kernel Learning algorithm
%
% This is the implementation of the algorithm described in
% "Infinite Kernel Learning" MPI Technical Report 178.
%
% Usage: 
%	[w, alpha, b, Theta, stat] = infkltrain(y, subproblem_func, C, opts)
%
% Input:
%	y = training labels [n x 1]
%	subproblem_func - function handle for the subproblem
%	C - regularing constant
%
% Optional Input:
%	opts - struct with the following possible fields (with
%	default values)
%		.maxiter (100) - Maximum number of iterations. One
%		iteration is the evaluation of the subproblem and
%		the masterproblem
%		.maxiter_ipopt (100) - Maximum number of iterations of
%		IPopt to solve the master problem	
%		.mkl_eps (1e-4) - precision for the MKL objective
%		.svm_eps (1e-8) - precision of the SVM solver
%		.decay_eps (1e-5) - weights below this value are
%		pruned away
%		.objective_eps (1e-5) - Convergence threshold. If
%		the relative (!) decrease of the objective value is
%		below this threshold, the program is stopped
%
% Output: 
%	w - mixing weights [M x 1]
%	[alpha,b] - SVM parameter of the solution
%	Theta - Kernel coefficients (as defined by subproblem function)
%	stat - struct with statistics of the algorithm (running time, termination reason, etc)
%
% Example: 
%	
%	see chessboard_example.m for an example run of this function
% 
%
% Peter Gehler 07/2008 pgehler@tuebingen.mpg.de

function [w, alpha, b, Theta, stat] = IKLmain(y, subproblem_func, C, opts)


%
% You need to update the path to IPOPT!
%
warning('off','MATLAB:dispatcher:pathWarning')
addpath('./libsvm-matlab');
switch computer 
 case 'GLNX86' 
  addpath('./ipopt/ipopt-3.3.5-mumps-i686/matlab');
case 'GLNXA64' % cluster
  addpath('./ipopt/ipopt-3.3.5-mumps-amd64/matlab');
end  
warning('on','MATLAB:dispatcher:pathWarning')


start_time = cputime;

if numel(C)~=1 | C<=0
    error('input arguement C must be positive scalar');
end

% sanity check
assert(max(y(:))==1);
assert(min(y(:))==-1);

maxiter_default = 100; % how many iterations ? 
maxiter_ipopt_default = 100; % how many (inner) iterations in Ipopt
mkl_eps_default = 1e-4; % precision for MKL solver (<svm_eps_default!!!)
svm_eps_default = 1e-6; % precision for SVM solver
decay_eps_default = 1e-5; % weights below this value will be removed
objective_eps_default = 1e-5; % threshold for the IKL algorithm

if ~exist('opts','var')
    opts = [];
end

if ~isfield(opts,'maxiter')
    opts.maxiter = maxiter_default;
end
if ~isfield(opts,'maxiter_ipopt')
    opts.maxiter_ipopt = maxiter_ipopt_default;
end
if ~isfield(opts,'mkl_eps')
    opts.mkl_eps = mkl_eps_default;
end
if ~isfield(opts,'svm_eps')
    opts.svm_eps = svm_eps_default;
end
if ~isfield(opts,'decay_eps')
    opts.decay_eps = decay_eps_default;
end
if ~isfield(opts,'objective_eps')
    opts.objective_eps = objective_eps_default;
end

DO_TEST = 1;


% needed for Ipopt
global infkl_f infkl_df infkl_w
global infkl_alpha infkl_b
global svm_counter_iteration deriv_counter

global infkl_multiplier % for sanity checks in the subfunction

% for statistics
svm_counter_iteration = 0;
svm_counter_total = 0;

old_infkl_f = Inf;
lbc = 1;
ubc = 1;

% obtain an initial starting point from the subproblem function
global Ks
[Ks,Theta,parm] = feval(subproblem_func,'init',[],[],[],[]);
%[Ks,Theta,parm] = feval(subproblem_func,[],[],[],[]);

nof_constraints = size(Ks,3);
nof_ipopt_died = 0;

w0 = ones(size(Ks,3),1)/size(Ks,3);

stats = [];
active_set = 1:numel(w0);

termination_reason = '';

for iter=1:opts.maxiter

    if any(isnan(Ks(:)))
	error('Gram matrix has NaN entry');
    end
    
    if (DO_TEST)
	old_alpha = infkl_alpha;
	old_b = infkl_b;
	old_w = w0(1:end-1);
    end

    nof_ipopt_tries = 0;
    max_ipopt_tries = 7;
    
    
    while nof_ipopt_tries < max_ipopt_tries
	nof_ipopt_tries = nof_ipopt_tries + 1;
	try
	    infkl_f = Inf;
	    infkl_df = Inf;
	    infkl_w = Inf(size(w0));

	    infkl_alpha = Inf;
	    infkl_b = Inf;

	    % 0 <= beta <= 1
	    lb = zeros(size(w0));
	    ub = Inf(size(w0));

	    svm_counter_iteration = 0;
	    deriv_counter = 0;

	    % call MKL with fininte number of kernels Ks
	    [w,multiplier,numiter] = ipopt(w0,lb,ub,lbc,ubc,...
					   @(w) compute_objective(w,y,Ks,C,opts.svm_eps),...
					   @(w) compute_gradient(w,y,Ks,C,opts.svm_eps),...
					   @(w) compute_constraints(w),...
					   @(w, rso) compute_jacobian(w, rso),...
					   '', [],'',[],...
					   'jac_c_constant','yes',...
					   'hessian_approximation','limited-memory',...
					   'mu_strategy','adaptive',...
					   'tol',opts.mkl_eps,...
					   'print_level',0,... % 5 for verbosity
					   'max_iter',opts.maxiter_ipopt);
	    break;
	catch
	    
	    assert(size(Ks,3)>1,'IpOpt died but #Gram matrices =1');

	    fprintf('Infkl: IpOpt Died, removing constraint and retrying\n');

	    %
	    % Ipopt died. Now we try the following heuristic:
            % remove one constraint and retrain. We remove the
            % constraint (gram matrix) which is closest in absolute
            % deviation from the next closest one
	    %
	    nof_ipopt_died = nof_ipopt_died + 1;
	    
	    diff = Inf(size(Ks,3),size(Ks,3));
	    for ii=1:(size(Ks,3)-1)
		for jj=(ii+1):size(Ks,3)
		    diff(ii,jj) = sum(sum(abs(Ks(:,:,ii)-Ks(:,:,jj))));
		end
	    end
	    assert(size(diff,1)==size(diff,2))
	    
	    [i1,i2] = ind2sub(size(diff),find(diff==min(diff(:))));

	    if w0(i1) > w0(i2)
		remove_ind = i2;
	    else
		remove_ind = i1;
	    end

	    % now we remove the constraint and restart IpOpt
	    if iscell(Theta)
		Theta{remove_ind} = [];
	    else
		Theta(remove_ind,:) = [];
	    end
	    Ks(:,:,remove_ind) = [];
	    w0(remove_ind) = [];
	    w0 = w0/sum(w0);
	    active_set(remove_ind) = [];

	    old_infkl_f = Inf;
	    
	end
    end
    if nof_ipopt_tries == max_ipopt_tries
	termination_reason = 'ipopt_died';
	error('IPOpt died too often');
	break;
    end
    
    numiter = numiter + 1; % C->Matlab
    %assert(numiter==deriv_counter,'#derivative calls != #Ipopt iterations');

    infkl_multiplier = multiplier;
    
    fprintf('IKL: IpOpt iter %d,',numiter);
    obj_decrease = (old_infkl_f-infkl_f)./infkl_f;
    fprintf('Calls F:%d, dF:%d,',svm_counter_iteration,deriv_counter);
    fprintf('obj: %.05g,',infkl_f);
    fprintf('rel.decrease: %.05f (thresh %.05f)\n',obj_decrease,opts.objective_eps);

    svm_counter_total = svm_counter_total + svm_counter_iteration;
    
    if (DO_TEST) 
	% it is possible that the objective function increases a
        % bit due to numerical errors. If two constraints are very
        % close IpOpt has a hard time to score the two kernel
        % matrices. Furthermore we are pruning out some weights
        % which might lead to precision error.
	assert(opts.maxiter_ipopt<=numiter || obj_decrease >= -1e-2,'Objective function increased');
    end
    old_infkl_f = infkl_f;

    % check if objective decrease is small
    if obj_decrease<opts.objective_eps
	termination_reason = 'obj_tol';
	break;
    end

    % add no new constraint, maximum of iterations reached
    if iter==opts.maxiter
	termination_reson = 'maxiter';
	break;
    end
    
    tmp_w = zeros(size(w,1),1);
    tmp_w(active_set,1) = w;
    
    %
    % now solve the SUBPROBLEM
    %
    %[newK,theta,parm] = feval(subproblem_func,infkl_alpha,multiplier.lambda,parm,tmp_w);
    [newK,theta,parm] = feval(subproblem_func,'search',infkl_alpha,multiplier.lambda,parm,tmp_w);
    

    if numel(newK)==0 % for some reason size([],3)==1, we need to
                   % catch this
	nof_new_constraints = 0;
    else
	nof_new_constraints = size(newK,3);
    end

    if iscell(theta)
	assert(numel(Theta)==nof_new_constraints);
    else
	assert(size(theta,1)==nof_new_constraints);
    end
    
    % ... found a new constraint? 
    if nof_new_constraints == 0
	termination_reason = 'no_constr';
	break;
    end
    if (DO_TEST)
	old_active_set = active_set;
    end

    % ... add them ...
    active_set = [active_set,nof_constraints + (1:nof_new_constraints)];
    nof_constraints = nof_constraints + nof_new_constraints;

    % ... and prune out decayed w ...
    decayed_w = find(w<opts.decay_eps);
    active_set(decayed_w) = []; % keep track of which points were selected
    
    % ... delete the corresponding Gram matrices ...
    Ks(:,:,decayed_w) = [];
    w(decayed_w) = [];
    Ks(:,:,end+(1:nof_new_constraints)) = newK;

    % ... restart at old point, add the last with 0 weight
    w0 = [w;zeros(nof_new_constraints,1)];
    
    % ... store parameters ...
    if iscell(theta)
	Theta{end+1} = theta;
	Theta{decayed_w} = [];
    else
	Theta(end+(1:nof_new_constraints),:) = theta;
	Theta(decayed_w,:) = [];
    end

    % ... remove precision error ... 
    % NO we can not do this, otherwise the objective function we
    % monitor might differ and could even increase solely due to
    % rounding errors. Therefore this is left commented out but the
    % assertion decrease>0 left in
    %w = w./sum(w);
    
    fprintf('current w: ');fprintf('%4f ',w); fprintf('\n');

    % ... and do assertions ...
    assert(size(Ks,3)==numel(active_set));
    assert(numel(active_set)==numel(unique(active_set)));
    assert(numel(w)+nof_new_constraints==numel(active_set));
    if iscell(Theta),assert(numel(Theta)==numel(w)+nof_new_constraints);
    else assert(size(Theta,1)==numel(w)+nof_new_constraints);end
    assert(size(Ks,3)==numel(w0));
    assert((opts.decay_eps~=Inf)|(all(active_set==1:numel(active_set))));
    
end


% report why algorithm terminated ...
fprintf('IKL: Algorithm terminated beracuse ''%s''\n',get_termination_reason(termination_reason));

% ... prune out decayed w ...
decayed_w = find(w<opts.decay_eps);
%active_set(decayed_w) = [];
w(decayed_w) = [];
assert(size(w,2)==1,'dimension error w should be [N,1]');

if iscell(Theta)
    Theta{decayed_w} = [];
    assert(numel(Theta)==numel(w),'#param != #mixture coeffs');
else
    Theta(decayed_w,:) = [];
    assert(size(Theta,1)==size(w,1),'#param != #mixture coeffs');
end

% ... remove precision error ... 
w = w./sum(w); % now we can do it.

% ... and print final statistics!
fprintf('final w = : ');fprintf('%.4f ',w);fprintf('\n');
fprintf('IKL: #subproblem calls: %d\n',iter);

% ... copy the final solution to output args
alpha = infkl_alpha; 
b = infkl_b;

stat.subproblem_calls = iter;
stat.termination_reason = get_termination_reason(termination_reason);
stat.svm_calls = svm_counter_total;
stat.nof_constraints = nof_constraints;
stat.total_time = start_time-cputime;
stat.nof_ipopt_died = nof_ipopt_died;

return;

% end of IKL function


%
% What follows now is a SimpleMKL implementation to solve the
% restricted master problem.
%

% objective function
function f = compute_objective(w,y,Ks,C,svm_eps)
global svm_counter_iteration 
svm_counter_iteration = svm_counter_iteration +1;
f = compute_func(w,y,Ks,C,svm_eps);
assert(numel(f)==1,'dimension mismatch');

% gradient function
function df = compute_gradient(w,y,Ks,C,svm_eps)
global deriv_counter
deriv_counter = deriv_counter + 1;
[foo,df] = compute_func(w,y,Ks,C,svm_eps);
assert(all(size(df)==size(w)),'df: dimension mismatch');

% constraints sum(w) = 1
function f = compute_constraints(w)
f = sum(w);

% Jacobian
function J = compute_jacobian(w,return_structure_only);
J = sparse(ones(size(w)))';


% compute the objective and gradient
% f = sum(alpha) - 0.5 sum( w_i * al' *K_i*al)
% df = -0.5 * al' * K * al
function [f,df] = compute_func(w,y,Ks,C,svm_eps)

global infkl_w infkl_f infkl_df infkl_alpha infkl_b
if any(size(infkl_w)~=size(w)) | any(infkl_w~=w)
    
    kernel.type = 'custom';
    [al,b] = libsvm(y,Kbeta(Ks,w,1),C,kernel,1,'svm_eps',svm_eps);

    infkl_alpha = al;
    infkl_b = b;
    
    svset = find(al~=0);

    assert(numel(svset)>0,'no support vectors selected');

    al = al(svset);
    for i=1:size(Ks,3)
	norm(i,1) = 0.5 * al'*Ks(svset,svset,i)*al;
    end
    infkl_f = sum(abs(al)) - sum(w.*norm);
    infkl_df = -norm;

    infkl_w = w;
end

f = infkl_f;
df = infkl_df;

if ~(numel(f)==1 || C>0 || numel(infkl_b)==1   ||     numel(infkl_alpha)==numel(y))
    fprintf('assert will fail\n');
    assert(0);
end

% sanity checks
assert(all(size(infkl_w)==size(infkl_df)));
assert(numel(f)==1);
assert(C>0);
assert(numel(infkl_alpha)==numel(y));
assert(numel(infkl_b)==1);


%
% Some termination reasons which could occur
%
function str = get_termination_reason(descr)
switch descr
 case 'obj_tol'
  str = 'Objective change smaller than threshold';
 case 'no_constr'
  str = 'No violating contraint found by subproblem';
 case 'maxiter'
  str = 'Maximum of iterations reached';
 case 'ipopt_died'
  str = 'IpOpt died unexpectedly';
 otherwise
  str = 'unkown reason of termination';
end

