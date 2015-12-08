function Hessian = HessAugMexbias(K,wj,activeset,v,C,cgamma,cgammab,Hessian,alpha1,alpha2,B)

len = length(activeset);
for j=1:len
	ind = activeset(j);
    Hessian = Hessian + alpha1(ind)*K(:,:,ind);
	Hessian = Hessian + alpha2(ind)*(v(:,ind)*v(:,ind)');
end;

if ~exist('B','var')
    Hessian = Hessian + cgammab;
else
    Hessian = Hessian + cgammab*B*B';
end;