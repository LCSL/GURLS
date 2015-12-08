% SimpleMKL implementation based on Ipopt and Libsvm
%
% Usage:
%	[w, alpha, b] = simplemkl(y, Ks, C, mkl_eps, svm_eps)
%
% Input:
%	y - vector of labels (n x 1)
%	Ks - M Gram matrices (n x n x M)
%	C - Regularizer
%
% optional Inputs:
%	mkl_eps - tolerance for MKL solver (default: 1e-5)
%	svm_eps - tolerance for SVM solver (default: 1e-6)
%		Note: it should be svm_eps<mkl_eps
%
% Outputs: 
%	w - weights for Gram matrix (M x 1)
%	alpha - Multipliers for SVM 
%	b - SVM threshold
%
% Example:
%	y = sign(randn(100,1));
%	for d=1:10 % 10 random Gram matrices
%	    X = randn(100,100);
%	    Ks(:,:,d) = X'*X;
%	end
%	[w,al,b] = simplemkl(y,Ks,100);
%
% You need to update the paths to the libsvm function and the ipopt solver 
%
% Peter Gehler 07/2008 pgehler@tuebingen.mpg.de

function [w, alpha, b,multiplier] = simplemkl(y, Ks, C, mkl_eps, svm_eps)

%
% You need to modify to your system!
%
warning('off','MATLAB:dispatcher:pathWarning')
addpath('./libsvm');
switch lower(computer)
 case 'glnx86' 
  addpath('./ipopt/ipopt-3.3.5-mumps-i686/matlab');
 case 'glnxa64' % cluster
  addpath('./ipopt/ipopt-3.3.5-mumps-amd64/matlab');
 otherwise

  error('unknown architecture');
end
warning('on','MATLAB:dispatcher:pathWarning')

% sanity checks
assert(size(Ks,1)==size(Ks,2),'no symmetric kern matrices');
assert(numel(y)==size(Ks,1),'dimension mismatch label gram matrix');
assert(C>0,'C<=0');
assert(max(y)>0,'no positive label');
assert(min(y)<0,'no negative label');
assert(~any(isnan(Ks(:))),'Gram matrix contains NaNs');
assert(~any(isinf(Ks(:))),'Gram matrix contains Infs');

if ndims(Ks)==2
    nGramMatrices = 1;
else
    nGramMatrices = size(Ks,3);
end

if nGramMatrices < 1
    error(['You must provide at least one kernel to mkltrain.']);
end

% initialize
w0 = ones(nGramMatrices,1)/nGramMatrices;

% needed for the Ipopt solver
global stored_norm stored_w stored_f
global objective_counter retrain_counter gradient_counter

if nargin < 4
    mkl_eps = 1e-5;
end
if nargin < 5
    svm_eps = 1e-6;
end

% 0 <= w <= Infty
lb = zeros(size(w0));
ub = Inf(size(w0));

% sum(w) = 1
lbc = 1;
ubc = 1;

% sometimes Ipopt dies, restarting with slightly different w's
% solves this problem and we do this at max 5 times
nof_max_tries = 5; 
nof_tries = 0;
while nof_tries<nof_max_tries
    try
	nof_tries = nof_tries+1;

	objective_counter = 0;
	retrain_counter = 0;
	gradient_counter = 0;

	stored_f = Inf;
	stored_w = Inf(size(w0));
	stored_norm = Inf;

	% SIMPLEMKL
	[w,multiplier] = ipopt(w0,lb,ub,lbc,ubc,...
		  @(w) compute_objective(w,y,Ks,C,svm_eps),...
		  @(w) compute_gradient(w,y,Ks,C,svm_eps),...
		  @(w) compute_constraints(w),...
		  @(w, rso) compute_jacobian(w, rso),...
		  '',[],'',[],...
		  'jac_c_constant','yes','hessian_approximation','limited-memory',...
		  'mu_strategy','adaptive',...
		  'print_level',0,... %5 for verbosity, 0 for silence
		  'tol',mkl_eps);
	break;
    catch
	% just add a bit of noise to w and restart
	w0 = w0 + 1e-6*randn(size(w0));
	fprintf('mkltrain_mysimple: IpOpt died: retrying\n');
    end
end

if nof_tries>=nof_max_tries
    % It did not converge!
    error('IpOpt died!');
end

% ... compute the final SVM parameters 
kernel.type = 'custom';
[alpha,b] = libsvm(y,Kbeta(Ks,w,1),C,kernel,1,'svm_eps',1e-5);

fprintf('MKL: calls F:%d,',objective_counter);
fprintf('SVM:%d,',retrain_counter);
fprintf('dF:%d\n',gradient_counter);

% end of function

%
% OBJECTIVE FUNCTION
%
function f = compute_objective(w,y,Ks,C,svm_eps)

global stored_w stored_f stored_norm
global objective_counter retrain_counter

objective_counter = objective_counter + 1;

if all(stored_w == w) && ~isinf(stored_f)
    f = stored_f;
    return;
end

retrain_counter = retrain_counter + 1;
kernel.type = 'custom';
indx = find(w~=0);
if ndims(Ks)==3
    [al,b] = libsvm(y,Kbeta(Ks(:,:,indx),w(indx),1),C,kernel,1,'svm_eps',svm_eps);
else
    [al,b] = libsvm(y,Kbeta(Ks,w,1),C,kernel,1,'svm_eps',svm_eps);
end    
ind = find(al~=0);
al = al(ind);
for i=1:numel(w)
    norm(i,1) = double(0.5 * al'*Ks(ind,ind,i)*al);
end
stored_norm = norm;
stored_w = w;

f = double(sum(abs(al)) - sum(w.*norm));
stored_f = f;

assert(~isnan(f),'objective is NaN');


%
% GRADIENT OF OBJECTIVE FUNCTION
%
function df = compute_gradient(w,y,Ks,C,svm_eps)

global stored_w stored_f stored_norm gradient_counter 
global retrain_counter;

gradient_counter = gradient_counter +1;


if all(stored_w == w)
    df = - stored_norm;
    return;
end

retrain_counter = retrain_counter + 1;
kernel.type = 'custom';
indx = find(w~=0);
if ndims(Ks)==3
    [al,b] = libsvm(y,Kbeta(Ks(:,:,indx),w(indx),1),C,kernel,1,'svm_eps',svm_eps);
else
    [al,b] = libsvm(y,Kbeta(Ks,w,1),C,kernel,1,'svm_eps',svm_eps);
end    
ind = find(al~=0);
al = al(ind);
for i=1:numel(w)
    norm(i,1) = double(0.5 * al'*Ks(ind,ind,i)*al);
end
stored_norm = norm;
stored_w = w;

f = double(sum(abs(al)) - sum(w.*norm));
stored_f = f;
 
df = -norm;

assert(~any(isnan(df)),'derivative is NaN');

%
% CONSTRAINT sum(w) =1
%
function f = compute_constraints(w)
f = sum(w);


%
% JACOBIAN
%
function J = compute_jacobian(w,return_structure_only);
J = sparse(ones(size(w)))';


