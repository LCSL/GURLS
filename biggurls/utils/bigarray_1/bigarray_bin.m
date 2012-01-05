classdef bigarray_bin < bigarray_file

%-----------------------------------------------------------------------------------------------------------------------

methods
function T = bigarray_bin(path, varargin)

T@bigarray_file(path, 'bin', varargin{:});

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Static, Hidden)
function T = loadobj(U)

T = U;

T.base = gl_abspath(T.path);

T.ext = 'bin';

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

if c == 0
    s = ReadHeader(fn);
    return;
end

[m, num] = T.BlockMap(c);

if isempty(field)
    s.num = num;
    s.var.(T.fields{1}) = m.Data.(T.fields{1});
else
    if ~strcmp(field, 'num'), error('invalid field'); end
    s = num;
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = protected)
function SaveBlock(T, c, s)

fn = T.FilePath(c);

if T.verb && (c > 0)
    fprintf('%s: writing\n', fn);
end

if c == 0
    path = fileparts(fn);
    if ~exist(path, 'file') && exist(fileparts(path), 'dir'), mkdir(path); end
    WriteHeader(fn, s);
    return;
end

h = fopen(fn, 'wb');
if h < 0, error('unable to open "%s"', fn); end

fwrite(h, uint64(s.num)      , 'uint64'    );
fwrite(h, s.var.(T.fields{1}), T.classes{1});

fclose(h);

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = protected)
function [m, num] = BlockMap(T, c)

fn = T.FilePath(c);

h = fopen(fn, 'rb');
if h < 0, error('unable to open "%s"', fn); end
num = double(fread(h, 1, 'uint64'));
fclose(h);

m = memmapfile(fn, 'Offset', 8, 'Format', {T.classes{1}, [T.sizes{1}, num], T.fields{1}}, 'Repeat', 1);

end
end

%-----------------------------------------------------------------------------------------------------------------------

end

%***********************************************************************************************************************

function s = ReadHeader(fn)

h = fopen(fn, 'r');
if h < 0, error('unable to open "%s"', fn); end

s = struct;

while true

    line = fgetl(h);
    if isequal(line, -1), break; end

    [field, str] = strtok(line, '=');
    str = str(2 : end);

    switch field
    case 'type'
        val = str;
    case 'blockSize'
        val = str2double(str);
    case {'fields', 'classes'}
        if isempty(str)
            val = {};
        else
            val = {str};
        end
    case 'sizes'
        if isempty(str)
            val = {};
        else
            val = textscan(str, '%f');
            val{1} = val{1}';
        end
    case 'auto'
        if isempty(str)
            val = [];
        else
            val = logical(str2double(str));
        end
    otherwise
        error('invalid field "%s"', field);
    end

    s.(field) = val;

end

fclose(h);

end

%***********************************************************************************************************************

function WriteHeader(fn, s)

h = fopen(fn, 'w');
if h < 0, error('unable to open "%s"', fn); end

fields = fieldnames(s);

for i = 1 : numel(fields)

    val = s.(fields{i});

    switch fields{i}
    case 'type'
        str = val;
    case 'blockSize'
        str = num2str(val);
    case {'fields', 'classes'}
        if isempty(val)
            str = '';
        else
            if numel(val) > 1, error('multiple fields not supported'); end
            str = val{1};
        end
    case 'sizes'
        if isempty(val)
            str = '';
        else
            if numel(val) > 1, error('multiple fields not supported'); end
            str = sprintf('%u ', val{1});
            str = str(1 : end - 1);
        end
    case 'auto'
        if isempty(val)
            str = '';
        else
            str = num2str(val);
        end
    otherwise
        error('invalid field "%s"', fields{i});
    end

    fprintf(h, '%s=%s\n', fields{i}, str);

end

fclose(h);

end