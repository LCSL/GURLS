function W = rp_projections(D,d,kernel)
switch kernel
    case 'gaussian'
        W = sqrt(2)*randn(d,D);
    otherwise
        error('cannot sample from that yet.')
end
end
