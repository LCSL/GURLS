function abs = gl_isabspath(path)

if isunix
    abs = ~isempty(strmatch('/', path));
else
    abs = ~isempty(regexp(path, '^[A-Za-z]:\', 'once'));
end

return;