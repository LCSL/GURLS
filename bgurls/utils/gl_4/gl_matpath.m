function p = gl_matpath(q)

persistent path;
if isempty(path), path = ''; end

p = path;

if nargin >= 1
    if ~ischar(q), error('invalid load path'); end
    if ~isempty(q), q = gl_abspath(q); end
    path = q;
end

return;