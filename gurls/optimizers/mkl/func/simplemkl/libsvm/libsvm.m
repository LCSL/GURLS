% [alpha,b,obj] = libsvm(y,K,C,kernel,l1)
%
% Interface to the LIBSVM for matlab
% In: 
%  'y' - labels
%  'K' - either Kernel matrix or training data
%  'C' - regularization parameter
% optional
%  'kernel' - struct with kernel defintion
%  'l1' - if 1 : l1 penalization, otherwise: l2 penalization
%
% out: 
%  'alpha','b' - svm parameter
%  'obj' - objective value
%
% prediction is sign(K*alpha+b)
%
% Peter Gehler, March 6th 2008
% pgehler@tuebingen.mpg.de
function [alpha,b,obj] = libsvm(y,K,C,kernel,l1,varargin)


if nargin < 3
    error('invalid number of arguments to libsvm');
end

%
% default values
%
if ~exist('l1','var')
    l1 = 1;
end

if size(y,1) == 1
    y = y';
end

if size(y,1) ~= size(K,1)
    error('Dimension mismatch between labels and Kern matrix');
end
if size(y,2) ~= 1
    error('invalid format of y');
end
if (size(K,1) ~= size(K,2)) || ndims(K)~=2
    error('invalid format of Kern matrix');
end

kernelType = 4; %Custom Kernel
cachesize = 40; % in MB
gamma = 0; % not needed for custom kernel
degree = 3; % degree for polynomial
svm_type = 0; %>1 for nu-svm
nu = 0; % nu-svm
p=0.001;
shrinking = 1;
balanced_ridge = 0;
coef0 = 1;
svm_eps = 1e-5;
weights = [];
if nargin > 5
    for i=1:2:nargin-5
	switch varargin{i}
	 case 'svm_eps',svm_eps = varargin{i+1};
	 case 'weights',weights = varargin{i+1};
	 otherwise 
	  error(['unknown argument ''',varargin{i},'''\n']);
	end
    end
end


if ~exist('kernel,','var') || ~numel(kernel)
    clear kernel
    kernel.type = 'custom';
end

switch kernel.type
 case 'rbf'
  gamma = kernel.gamma;
 case 'poly' % 1 -- polynomial: (gamma*u'*v + coef0)^degree
  degree = kernel.degree;
  gamma = kernel.gamma;
  coef0 = kernel.coef0;
  degree = kernel.degree;
 case 'custom'
  x = repmat( (1:size(K,1))',1,2);
 otherwise
  error('unknown kernel');
end

number_classes = numel(unique(y));

if number_classes > 2 % Multiclass : 1 vs Rest

    classes = unique(y);
    for cl = 1:number_classes
	trY = y;
	trY( y~=classes(cl) ) = -1;
	trY( y==classes(cl) ) = 1;

	assert( sum((trY==1)|(trY==-1)) == numel(trY));
	
	[alpha(:,cl),b(cl)] = libsvm_train(x,trY,C,svm_eps,kernel,gamma,p,shrinking,degree,coef0,K,l1,weights);
	
	if nargout > 2
	    obj(cl,1) = libsvm_objective(trY,K,alpha(:,cl),b(cl),C,l1);
	end

    end

    assert(size(alpha,2) == number_classes);
    assert(numel(b) == number_classes);
    if nargout > 2
	assert(numel(obj)==number_classes);
    end
    
elseif number_classes == 2 % binary
    
    if min(y(:)) ~= -1 || max(y(:))~=1
	error('labels must be in set {-1,+1}');
    end
    
    [alpha,b] = libsvm_train(x,y,C,svm_eps,kernel,gamma,p,shrinking,degree,coef0,K,l1,weights);

    assert(numel(alpha) == size(K,1));
    assert(numel(b) == 1);
    
    if nargout > 2
	obj = libsvm_objective(y,K,alpha,b,C,l1);
    end

else
    error('invalid number of classes in y');
end




function [alpha,b] = libsvm_train(x,y,C,svm_eps,kernel,gamma,p,shrinking,degree,coef0,K,l1,weights)

assert(size(x,1)==size(y,1));
assert(size(y,2)==1);
assert(size(K,1)==size(K,2));
assert(size(K,1)==numel(y));
assert(isfield(kernel,'type')>0);

switch kernel.type
 case 'rbf'
  kernelType = 2;
 case 'poly' % 1 -- polynomial: (gamma*u'*v + coef0)^degree
  kernelType = 1;
 case 'custom'
  kernelType = 4;
 otherwise
  error('unknown kernel');
end



if numel(weights)==0
    switch kernel.type
     case 'custom'
      if l1
	  [al,xSV,bias0]=libsvm_classifier({'X',double(x)},{'Y',y}, ...
					   {'C',C},...
					   {'eps',svm_eps},{'kerneltype',kernelType},{'gamma',gamma},...
					   {'p',p},{'shrinking',shrinking},... 
					   {'degree',degree},...
					   {'coef0', coef0},...
					   {'kmatrix',double(K)});
      else
	  [al,xSV,bias0]=libsvm_classifier({'X',double(x)},{'Y',y}, ...
					   {'C',1e6},...
					   {'eps',svm_eps},{'kerneltype',kernelType},{'gamma',gamma},...
					   {'p',p},{'shrinking',shrinking},... 
					   {'degree',degree},...
					   {'coef0', coef0},...
					   {'kmatrix',double(K)+eye(size(K,1))/C});
      end
     otherwise
      if l1
	  [al,xSV,bias0]=libsvm_classifier({'X',double(K)},{'Y',y}, ...
					   {'C',C},...
					   {'eps',svm_eps},{'kerneltype',kernelType},{'gamma',gamma},...
					   {'p',p},{'shrinking',shrinking},... 
					   {'degree',degree},...
					   {'coef0', coef0});
      else
	  error('currently not supported');
      end
      
    end
else
    fprintf('using weights!\n');
    
    switch kernel.type
     case 'custom'
      if l1
	  [al,xSV,bias0]=libsvm_classifier_weight({'X',double(x)},{'Y',y}, ...
						  {'C',C},...
						  {'W',weights},...
						  {'eps',svm_eps},{'kerneltype',kernelType},{'gamma',gamma},...
						  {'p',p},{'shrinking',shrinking},... 
						  {'degree',degree},...
						  {'coef0', coef0},...
						  {'kmatrix',double(K)});
      else
	  [al,xSV,bias0]=libsvm_classifier_weight({'X',double(x)},{'Y',y}, ...
						  {'C',1e6},...
						  {'W',weights},...
						  {'eps',svm_eps},{'kerneltype',kernelType},{'gamma',gamma},...
						  {'p',p},{'shrinking',shrinking},... 
						  {'degree',degree},...
						  {'coef0', coef0},...
						  {'kmatrix',double(K)+eye(size(K,1))/C});
      end
     otherwise
      if l1
	  [al,xSV,bias0]=libsvm_classifier_weight({'X',double(K)},{'Y',y}, ...
						  {'C',C},...
						  {'W',weights},...
						  {'eps',svm_eps},{'kerneltype',kernelType},{'gamma',gamma},...
						  {'p',p},{'shrinking',shrinking},... 
						  {'degree',degree},...
						  {'coef0', coef0});
      else
	  error('currently not supported');
      end
      
    end


end

assert(max(xSV(:)) <= numel(y));

b = bias0 * y(1); 
al = al * y(1);
alpha = zeros(numel(y),1);
alpha(xSV(:,1)) = al;

function [obj] = libsvm_objective(y,K,alpha,b,C,l1);

if numel(K) == 0
%if ~strcmp(kernel.type,'custom')
    warning('objective calculation only supported for custom kernel');
    obj = Inf;
    return;
end

Kalpha = K*alpha;
out = Kalpha + b;
sv = find(out<1);

if l1
    obj = 0.5 * (sum(abs(alpha))/C + sum( abs((1-y(sv).*out(sv)))));
else
    obj = 0.5 * (alpha'*Kalpha/C + sum( (1-y(sv).*out(sv)).^2));
end

