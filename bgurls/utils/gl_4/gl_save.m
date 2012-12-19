function gl_save(path, varargin)

path = gl_abspath(path);

if (numel(varargin) == 1) && isstruct(varargin{1})
    s = varargin{1}; %#ok
elseif isempty(varargin)
    list = evalin('caller', 'whos');
    list = {list(~[list.global]).name};
    s = struct;
    for i = 1 : numel(list)
        s.(list{i}) = evalin('caller', list{i});
    end
elseif iscellstr(varargin)
    s = struct;
    for i = 1 : numel(varargin)
        s.(varargin{i}) = evalin('caller', varargin{i});
    end
else
    error('invalid arguments');
end

current = gl_matpath;
try
    gl_matpath(fileparts(path));
    save(path, '-struct', 's');
    err = [];
catch err
end
gl_matpath(current);

if ~isempty(err), rethrow(err); end

return;