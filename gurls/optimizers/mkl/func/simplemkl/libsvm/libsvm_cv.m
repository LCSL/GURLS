% [alpha,b,obj,optimalC] = libsvm_cv(y,K,Cs,nFolds,kernel,l1)
%
% Cross Validtion Interface to the LIBSVM for matlab
% In: 
%  'y' - labels
%  'K' - either Kernel matrix or training data
%  'Cs' - list of regularization parameters
%  'nFolds' - number of Cross Validation splits to do
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
function [alpha,b,obj,optimalC] = libsvm_cv(y,K,Cs,nFolds,kernel,l1,varargin)


%
% default values
%

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

if numel(unique(y))~=2
    error('must give exactly two classes');
end


[tr_split,te_split] = create_cvsplit(y,nFolds);

if numel(Cs) > 1
    for cind = 1:numel(Cs)
	for f=1:nFolds
	    if exist('l1','var')
		[alpha,b] = libsvm(y(tr_split{f}),K(tr_split{f},tr_split{f}),Cs(cind),kernel,l1);
	    elseif exist('kernel','var')
		[alpha,b] = libsvm(y(tr_split{f}),K(tr_split{f},tr_split{f}),Cs(cind),kernel);
	    else
		[alpha,b] = libsvm(y(tr_split{f}),K(tr_split{f},tr_split{f}),Cs(cind));
	    end
	    
	    Ypred = sign(K(te_split{f},tr_split{f})*alpha + b);
	    
	    err(cind,f) = mean(Ypred ~= y(te_split{f}));
	end
    end
    
    err = err + 1e-9 * randn(size(err)); % break ties
    
    err = mean(err,2);
    assert(numel(err) == numel(Cs));
    
    [minErr,bestInd] = min(err);
    
    optimalC = Cs(bestInd);

else
    optimalC = Cs;
end


if exist('l1','var')
    [alpha,b,obj] = libsvm(y,K,optimalC,kernel,l1);
elseif exist('kernel','var')
    [alpha,b,obj] = libsvm(y,K,optimalC,kernel);
else
    [alpha,b,obj] = libsvm(y,K,optimalC);
end

