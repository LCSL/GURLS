function [opt] = set_sigma(opt, sigma)
% insert sigma to the slot 'opt.paramsel.sigma'
if ~isprop(opt,'paramsel')
    opt.newprop('paramsel', struct());
    % signature
    opt.paramsel.manual_sigma = true;
end
opt.paramsel.sigma = sigma;
end
