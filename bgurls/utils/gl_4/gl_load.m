function r = gl_load(path, varargin)

path = gl_abspath(path);

current = gl_matpath;
try
    gl_matpath(fileparts(path));
    s = load(path, varargin{:});
    err = [];
catch err
end
gl_matpath(current);

if ~isempty(err), rethrow(err); end

if nargout == 0
    fields = fieldnames(s);
    for i = 1 : numel(fields)
        assignin('caller', fields{i}, s.(fields{i}));
    end
else
    r = s;
end

return;