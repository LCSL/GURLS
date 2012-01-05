classdef bigarray_mat < bigarray_file

%-----------------------------------------------------------------------------------------------------------------------

methods
function T = bigarray_mat(path, varargin)

T@bigarray_file(path, 'mat', varargin{:});

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Static, Hidden)
function T = loadobj(U)

T = U;

T.base = gl_abspath(T.path);

T.ext = 'mat';

T.Construct;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Hidden)
function display(T)

T.Display;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = protected)
function s = LoadBlock(T, c, field)

if nargin < 3, field = ''; end

fn = T.FilePath(c);

if T.verb && (c > 0) && isempty(field)
    fprintf('%s: reading\n', fn);
end

if isempty(field)
    s = load(fn);
else
    s = load(fn, field);
    s = s.(field);
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = protected)
function SaveBlock(T, c, s) %#ok

fn = T.FilePath(c);

if T.verb && (c > 0)
    fprintf('%s: writing\n', fn);
end

if c == 0
    path = fileparts(fn);
    if ~exist(path, 'file') && exist(fileparts(path), 'dir'), mkdir(path); end
    save(fn, '-struct', 's');
else
    save(fn, '-struct', 's', '-v7.3');
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

end