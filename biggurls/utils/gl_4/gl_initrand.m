function seed = gl_initrand(seed)

if nargout == 0
    if nargin < 1
        now = clock;
        seed = etime(now, [now(1) 1 1 0 0 0]) * 50;
    end
    if isstruct(seed)
        rand ('state', seed.u);
        randn('state', seed.n);
    else
        rand ('state', seed);
        randn('state', seed);
    end
else
    if nargin == 1, error('must either set or get seed'); end
    seed = struct;
    seed.u = rand ('state');
    seed.n = randn('state');
end

return;