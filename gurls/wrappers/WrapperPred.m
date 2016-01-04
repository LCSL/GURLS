function [y, opt] = WrapperPred(X, opt)
%
% WrapperPred(X, opt) gives the prediction for data in X using the fields
% in struct OPT.
%
% Use: For evaluating IterativeRLSWrapper functions.
%
% NOTE: This is function is probably obsolete and should be removed.

if ~isa(opt, 'GurlsOptions')
    % warning('Compatibility mode with GURLS 1.0');
    opt = GurlsOptions(opt);
end
if ~isprop(opt, 'rls')
    opt.newprop('rls', struct());
end

if isprop(opt, 'C')
    % kernelfun is an outdated option: should be removed from opt
    % opt.newprop('kernel', opt.kernelfun(opt.X, [], opt));
    
    % overwritting (any) values in .rls by those under main opt fields
    opt.rls.X = opt.X;
    opt.rls.C = opt.C;
    if strcmp(opt.kernel.type, 'linear')
        opt.rls.W = opt.rls.X'*opt.rls.C;
    end
    
    opt.newprop('predkernel', predkernel_traintest(X, [], opt));
    
    y = pred_dual(X, [], opt);
    
elseif isprop(opt, 'W')
    opt.rls.W = opt.W;
    y = pred_primal(X, [], opt);
else
    error('option structure does not contain estimator')
end