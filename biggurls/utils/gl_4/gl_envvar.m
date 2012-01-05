function var = gl_envvar(name, dflt)

if ~isempty(which(name))
    var = feval(name);
elseif nargin == 2
    var = dflt;
else
    error('function "%s" not found', name);
end    

return;