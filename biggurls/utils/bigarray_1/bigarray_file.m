classdef bigarray_file < bigarray

%-----------------------------------------------------------------------------------------------------------------------

properties (Access = protected)
    path;
end

properties (Transient, Access = protected)
    base;
    ext;
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Static)
function T = Obj(path, ext, varargin)

if nargin < 2, ext = ''; end

if isempty(ext)
    [ans, ans, ext] = FindBlocks2(gl_abspath(path), '');
    if isempty(ext), ext = 'mat'; end
end

switch ext
case 'mat', T = bigarray_mat(path, varargin{:});
case 'bin', T = bigarray_bin(path, varargin{:});
otherwise , error('invalid type');
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = protected)
function T = bigarray_file(path, ext, varargin)

T@bigarray(varargin{:});

[T.base, T.path] = gl_abspath(path);

T.ext = ext;

T.Construct;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Hidden)
function U = saveobj(T)

U = T.saveobj@bigarray;

U.path = gl_savepath(U.base, U.path);

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function path = Path(T)

path = T.path;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function base = Base(T)

base = T.base;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = protected)
function Display(T)

T.Display@bigarray(T.path);

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function Move(T, newPath, over)

if nargin < 3, over = false; end

if T.read, error('read only'); end

T.Flush;

[base, path] = gl_abspath(newPath);
if strcmp(base, T.base), error('source and destination are the same'); end

b = bigarray.Obj(newPath, T.ext);
if (b.State > 0) && ~over, error('destination files already exist, overwrite not specified'); end
b.Clear;
delete(b);

list = T.FindBlocks;
for i = 1 : numel(list)
    fn1 = list{i};
    fn2 = FilePath2(base, T.ext, T.digs, i - 1);
    movefile(fn1, fn2);
end

T.path = path;
T.base = base;

T.Construct;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = protected)
function [list, digs] = FindBlocks(T)

[list, digs] = FindBlocks2(T.base, T.ext);

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = protected)
function DeleteBlocks(T, head)

list = T.FindBlocks;

for i = 1 : numel(list)
    if i == 1
        if head, delete(list{i}); end
    else
        if T.verb, fprintf('%s: deleting\n', list{i}); end
        delete(list{i});
    end
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = protected)
function fn = FilePath(T, c)

fn = FilePath2(T.base, T.ext, T.digs, c);

end
end

%-----------------------------------------------------------------------------------------------------------------------

end

%***********************************************************************************************************************

function [list, digs, ext] = FindBlocks2(base, ext)

[path, name, name2] = fileparts(base);
name = [name name2];

list = dir(sprintf('%s.*', base));
digs = [];

inds = false(1, numel(list));

for i = 1 : numel(list)

    fn = list(i).name(numel(name) + 2 : end);
    if sum(fn == '.') ~= 1, continue; end
    if (fn(1) == '.') || (fn(end) == '.'), error('invalid files found'); end
    [num, rest] = strtok(fn, '.');
    rest = rest(2 : end);
    if ~all(isstrprop(num, 'digit')), continue; end
    if list(i).isdir, error('invalid files found'); end
    inds(i) = true;

    if isempty(digs)
        digs = numel(num);
    elseif numel(num) ~= digs
        error('invalid files found');
    end

    if isempty(ext)
        ext = rest;
    elseif ~strcmp(rest, ext)
        error('invalid files found');
    end

end

list = sort({list(inds).name});

for i = 1 : numel(list)
    list{i} = fullfile(path, list{i});
    if ~strcmp(list{i}, FilePath2(base, ext, digs, i - 1)), error('invalid files found'); end
end

end

%***********************************************************************************************************************

function fn = FilePath2(base, ext, digs, c)

fn = sprintf('%s.%0*u.%s', base, digs, c, ext);

end