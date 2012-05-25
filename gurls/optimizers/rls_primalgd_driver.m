function gd = rls_primalgd_driver(X,y,opt)

XtX = X'*X; % d x d matrix.
Xty = X'*y; % d x T matrix.

gd = opt.gd;
if (opt.gd.method == 1)
	alpha1 = opt.gd.alpha1;
	nu = opt.gd.nu;
end
i = opt.gd.iter;
W = opt.gd.W;

if (opt.gd.method == 0)
    W = W + opt.paramsel.eta*(Xty - XtX*W);
elseif (opt.gd.method == 1)
  u=((i-1)*(2*i-3)*(2*i+2*nu-1))/((i+2*nu-1)*(2*i+4*nu-1)*(2*i+2*nu-3));
  w=4*(((2*i+2*nu-1)*(i+nu-1)) /((i+2*nu-1)*(2*i+4*nu-1)) );
  alpha2 = alpha1;
  alpha1 = W;
  W = alpha1 + u*(alpha1 - alpha2) +(w*opt.paramsel.eta)*(Xty - XtX*alpha1);
else
  error('invalid opt.gd.method')
end
gd.W = W;
if (opt.gd.method == 1)
	gd.alpha1 = alpha1;
end
