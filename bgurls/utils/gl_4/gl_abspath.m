function [use, store] = gl_abspath(in, current)

% TODO: 'store' output should not have (env vars / ~) expanded?
% TODO: 2nd argument '%check'?

if ischar(in)

    path = Replace(in);
    if isempty(path), error('empty path'); end

    if ~gl_isabspath(path)
        if (nargin < 2) || isempty(current)
            current = cd;
        elseif strcmp(current, '%none')
            error('"%s" is not an absolute path', in);
        end
        path = Concat(current, path);
    end

    use   = path;
    store = path;

elseif iscellstr(in) && ~isempty(in) && (numel(in) <= 2)

    if numel(in) == 1, in = {'.', in{1}}; end

    path = Replace(in{2});
    if isempty(path), error('empty path'); end

    if gl_isabspath(path)

        use   = path;
        store = path;

    elseif isempty(in{1})

        base = gl_matpath;
        if isempty(base)
            error('"%s": base path not provided (if loading/saving, use gl_load/gl_save)', path);
        end

        use   = Concat(base, path);
        store = {'', path};

    else

        base = Replace(in{1});

        if ~gl_isabspath(base)
            if (nargin < 2) || isempty(current)
                current = cd;
            elseif strcmp(current, '%none')
                error('"%s" is not an absolute path', in{1});
            end
            base = Concat(current, base);
        end

        use   = Concat(base, path);
        store = {'', path};

    end
    
else

    error('invalid path');

end

return;

%***********************************************************************************************************************

function path = Replace(path)

if isunix
    path = strrep(path, '\', '/');
else
    path = strrep(path, '/', '\');
end

if isunix && (strcmp(path, '~') || ~isempty(strmatch('~/', path)))
    path = path(2 : end);
    if ~isempty(path), path = path(2 : end); end
    if isempty(path)
        path = getenv('HOME');
    else
        path = fullfile(getenv('HOME'), path);
    end
    return;
end

[i, j] = regexp(path, '^%\w+', 'once');
if isempty(i), return; end

var  = path(2 : j);
rest = path(j + 1 : end);
if ~isempty(rest) && (rest(1) == filesep), rest = rest(2 : end); end

path = gl_envvar(var);
if ~gl_isabspath(path), error('"%s" is not an absolute path', path); end

if ~isempty(rest), path = fullfile(path, rest); end

return;

%***********************************************************************************************************************

function path = Concat(base, path)

while ~isempty(path) && (path(1) == '.')
    if strcmp(path, '.') || ~isempty(strmatch(['.' filesep], path))
        path = path(2 : end);
        if ~isempty(path), path = path(2 : end); end
    elseif strcmp(path, '..') || ~isempty(strmatch(['..' filesep], path))
        path = path(3 : end);
        if ~isempty(path), path = path(2 : end); end
        base = fileparts(base);
    else
        error('"%s": invalid path', path);
    end
end

if isempty(path)
    path = base;
else
    path = fullfile(base, path);
end

return;