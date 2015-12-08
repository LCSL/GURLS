%
% Test method for the subproblem function, compares C and Matlab code
%
function subproblem_product_objective_test

nFac_list = [1,10];
nPoints_list = [1,100,500];
sig_list = [0.1,1,5];


Ds = randn(50,50,5);
al = randn(50,1);
fac = randn(5,1);

try
    [f,df,H] = subproblem_product_objective(Ds,fac,al,0.1,0.1);
    assert(false,'wrong number of input arguments');
end

try
    [f,df,H] = subproblem_product_objective(Ds,fac);
    assert(false,'wrong number of input arguments');
end

try
    [f,df,H] = subproblem_product_objective(Ds);
    assert(false,'wrong number of input arguments');
end

try
    [f,df,H,foo] = subproblem_product_objective(Ds,fac,al);
    assert(false,'wrong number of output arguments');
end

try
    subproblem_product_objective(Ds,fac,al);
    assert(false,'wrong number of output arguments');
end

try
    [f,df,H] = subproblem_product_objective(Ds,fac,al,-sig);
    assert(false,'negative variance should fail');
end




for ii=1:numel(nFac_list)
    nFacs = nFac_list(ii);
    for jj=1:numel(nPoints_list)
	nPoints = nPoints_list(jj);

	clear D Ds fac al;
	
	for k=1:nFacs
	    D = rand(nPoints,nPoints);
	    D = D+D';
	    D = D-diag(diag(D));
	    assert(sum(diag(D))==0);
	    assert(sum(sum(abs(D-D')))==0);
	    
	    Ds(:,:,k) = D;
	end

	al = randn(nPoints,1);
	fac = rand(nFacs,1);

	

	[f,df,H] = subproblem_product_objective(Ds,fac,al);
	[f2,df2] = subproblem_product_objective(Ds,fac,al);
	[f3] = subproblem_product_objective(Ds,fac,al);
	
	[f4,df4,H4] = subprob_obj(Ds,fac,al);


	
	assert(f==f2);
	assert(f==f3);
	assert(abs(f-f4)<1e-5);
	assert(all(df==df2));
	assert(sum(abs(df-df4))<1e-5);
	assert(sum(sum(abs(H-H4)))<1e-5);
	

	% smoothed variant
	
	for kk=1:numel(sig_list)
	    sig = sig_list(kk);

	    
	    [f,df,H] = subproblem_product_objective(Ds,fac,al,sig);
	    [f2,df2] = subproblem_product_objective(Ds,fac,al,sig);
	    [f3] = subproblem_product_objective(Ds,fac,al,sig);
	    [f4,df4,H4] = subprob_obj_sig(Ds,fac,al,sig);

	    assert(f==f2);
	    assert(f==f3);
	    assert(abs((f-f4)/f)<1e-4);
	    assert(all(df==df2));
	    if sum(df)==0
		assert(sum(df2)==0);
	    else
		assert(sum(abs(df-df2))/sum(abs(df))<1e-4);
	    end
	    if sum(H(:))==0
		assert(sum(H4(:))==0);
	    else
		assert(sum(sum(abs((H-H4))))/sum(abs(H(:)))<1e-4);
	    end
	    
	    
	    [f,df,H] = subprob_obj(Ds,fac,al);
	    [f2,df2,H2] = subprob_obj_sig(Ds,fac,al,0);

	    assert(f==f2);
	    
	    
	end	
	
    end 
end

fprintf('Test passed\n');

function [f,df,H] = subprob_obj(Ds,fac,al)

D = 0;
for k=1:size(Ds,3)
    D = D + fac(k) * Ds(:,:,k);
end
K = exp(-D);
f = - 0.5*al'*K*al;

for k=1:numel(fac)
    tmp = Ds(:,:,k) .* K;
    df(k,1) = 0.5 * al' * tmp *al;
    for k2=1:k
	H(k,k2) = -0.5*al'*(tmp.*Ds(:,:,k2))*al;
    end
    
end


function [f,df,H] = subprob_obj_sig(Ds,fac,al,sig)

assert(sig>=0);

if sig==0
    scale  = 1;
else
    scale = (pi*sig).^(-(size(Ds,3)+1)/2);
end


D = [];
for k=1:size(Ds,3)
    D(:,:,k) = fac(k) * (Ds(:,:,k)./(1+sig*Ds(:,:,k)));
end
K = exp(-sum(D,3));
s = (prod(sqrt(1+sig*Ds),3));
tmp = K./s;

f = - 0.5*al'*tmp*al;

for k=1:numel(fac)
    tmp = Ds(:,:,k) ./ s .* K;
    df(k,1) = 0.5 * al' * tmp *al;
    for k2=1:k
	H(k,k2) = -0.5*al'*(tmp.*Ds(:,:,k2))*al;
    end
    
end

if (0)
f = f*scale;
df = df*scale;
H = H*scale;
end
