classdef bigarray < handle

%-----------------------------------------------------------------------------------------------------------------------

properties (Access = private)

    opts;

    perm;

end

properties (Transient, Access = protected)

    read;
    verb;

    digs;
    blockSize;
    blockSize2;

    fields;
    classes;
    sizes;
    auto;

end

properties (Transient, Access = private)

    parts; % number of blocks
    c;     % index of currently loaded block
    mod;   % has the currently loaded block been modified?
    cnt;   % number of items in the currently loaded block
    buf;   % currently loaded block

    n;     % number of items

end

%-----------------------------------------------------------------------------------------------------------------------

methods (Static)
function T = Obj(path, ext, varargin)

if nargin < 2, ext = ''; end

if isempty(path)
    if ~ismember(ext, {'', 'mem'}), error('missing path'); end
    T = bigarray_mem(varargin{:});
else
    T = bigarray_file.Obj(path, ext, varargin{:});
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = protected)
function T = bigarray(opts)

if nargin < 1, opts = ''; end    

T.opts = opts;

T.perm = [];

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function Refresh(T)

T.Construct;

% Note we do not call T.Load(0).  Any unwritten changes are discarded.

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = protected)
function Construct(T)

T.SetFlags;

[list, T.digs] = T.FindBlocks;

if isempty(list)

    T.blockSize = [];
    T.fields    = {};
    T.classes   = {};
    T.sizes     = {};
    T.auto      = [];

else

    s = T.LoadBlock(0);
    if isfield(s, 'type')
        % Type check added 2011/05/12.  It's inside a conditional so we
        % don't invalidate existing saved bigarrays.
        if ~strcmp(s.type, class(T)), error('type mismatch'); end
    end

    T.blockSize = s.blockSize;
    T.fields    = s.fields;
    T.classes   = s.classes;
    T.sizes     = s.sizes;
    T.auto      = s.auto;

end

T.SetBlockSize;

T.parts = max(0, numel(list) - 1);
T.c     = 0;
T.mod   = false;
T.cnt   = inf;
T.buf   = [];

if T.parts == 0
    T.n = 0;
else
    num = T.LoadBlock(T.parts, 'num');
    T.n = (T.parts - 1) * T.blockSize2 + num;
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = private)
function SetFlags(T)

T.read = any(T.opts == 'r');
T.verb = any(T.opts == 'v');

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = private)
function SetBlockSize(T)

if T.State == 0

    T.blockSize2 = [];

elseif T.blockSize > 0

    T.blockSize2 = T.blockSize;

elseif T.State < 2

    T.blockSize2 = [];

else

    siz = 0;

    for i = 1 : numel(T.fields)

        bytes = Bytes(T.classes{i});
        if bytes == 0, error('class "%s" requires explicit blockSize', T.classes{i}); end

        siz = siz + prod(T.sizes{i}) * bytes;

    end

    T.blockSize2 = ceil(-T.blockSize / siz);

end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function prev = Opts(T, opts)

prev = T.opts;

if nargin < 2, return; end

T.Load(0);

if isempty(opts)
    T.opts = '';
else
    switch opts(1)
    case '+' , T.opts = unique([T.opts, opts(2 : end)]);
    case '-' , T.opts = T.opts(~ismember(T.opts, opts(2 : end)));
    otherwise, T.opts = unique(opts);
    end
end

T.SetFlags;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function Init(T, blockSize, digs)

if T.read, error('read only'); end
if T.State > 1, error('fields already exist'); end

if (nargin < 2) || isempty(blockSize)
    blockSize = -50 * 1024 * 1024;
elseif blockSize == 0
    error('invalid blockSize');
end

if (nargin < 3) || isempty(digs)
    digs = 8;
elseif digs == 0
    error('invalid number of digits');
elseif digs < 0
    digs = numel(num2str(-digs));
end

T.digs      = digs;
T.blockSize = blockSize;

T.SetBlockSize;

T.DeleteBlocks(true);

T.SaveHeader;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function Append(T, d, varargin)

if T.read, error('read only'); end

if T.State == 0
    T.Init;
end

if isa(d, 'bigarray')
    T.AppendBig(d, varargin{:});
else
    T.AppendData(d, varargin{:});
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = private)
function AppendBig(T, U, varargin)

if ~isempty(varargin) && (isnumeric(varargin{1}) || islogical(varargin{1}))
    inds = varargin{1};
    if islogical(inds)
        inds = find(inds);
    elseif ~isempty(inds)
        if any(diff(inds) < 1), error('indices must be ascending'); end
    end
    if ~isempty(inds)
        if (inds(1) < 1) || (inds(end) > U.NumItems), error('invalid indices'); end
    end
    varargin(1) = [];
else
    inds = 1 : U.NumItems;
end

if ~isempty(varargin) && isa(varargin{1}, 'function_handle')
    f = varargin{1};
    varargin(1) = [];
else
    f = @Null;
end

if isempty(inds)
    T.CheckData(f(U.Read([], varargin{:})));
    return;
end

j2 = 0;
while j2 < numel(inds)

    j1 = j2 + 1;
    if T.cnt < T.blockSize2
        num = min(numel(inds) - j1 + 1, T.blockSize2 - T.cnt);
    else
        num = min(numel(inds) - j1 + 1, T.blockSize2);
    end
    j2 = j1 + num - 1;

    T.AppendData(f(U.Read(inds(j1 : j2), varargin{:})));

end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = private)
function AppendData(T, d)

[d, cnt] = T.CheckData(d);
if cnt == 0, return; end

T.Load(T.parts);

j2 = 0;
while j2 < cnt

    if T.cnt >= T.blockSize2
        T.Load(T.parts + 1);
        T.parts = T.parts + 1;
    end

    j1 = j2 + 1;
    num = min(T.blockSize2 - T.cnt, cnt - j1 + 1);
    j2 = j1 + num - 1;

    for i = 1 : numel(T.fields)
        colons = T.Colons(i);
        T.buf.(T.fields{i})(colons{:}, T.cnt + 1 : T.cnt + num) = d.(T.fields{i})(colons{:}, j1 : j2);
    end

    T.cnt = T.cnt + num;
    T.mod = true;

end

T.n = T.n + cnt;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = private)
function [d, cnt] = CheckData(T, d)

auto = ~isstruct(d);
if auto
    temp = struct;
    temp.v = d;
    d = temp;
end

fields = fieldnames(d);

if T.State < 2
    if isempty(fields), error('no fields'); end
    T.fields  = fields;
    T.classes = cell(1, numel(fields));
    T.sizes   = cell(1, numel(fields));
    T.auto    = auto;
    for i = 1 : numel(fields)
        siz = size(d.(fields{i}));
        if i == 1
            cnt = siz(end);
        elseif siz(end) ~= cnt
            error('size mismatch');
        end
        T.classes{i} = class(d.(fields{i}));
        T.sizes  {i} = siz(1 : end - 1);
    end
    T.SetBlockSize;
    T.SaveHeader;
else
    if numel(fields) ~= numel(T.fields), error('field mismatch'); end
    for i = 1 : numel(T.fields)
        siz = size(d.(T.fields{i}));
        siz(end + 1 : numel(T.sizes{i}) + 1) = 1;
        if i == 1
            cnt = siz(end);
        elseif siz(end) ~= cnt
            error('size mismatch');
        end
        if ~strcmp(class(d.(T.fields{i})), T.classes{i}), error('class mismatch'); end
        if ~isequal(siz(1 : end - 1), T.sizes{i}), error('size mismatch'); end
    end
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function Flush(T)

T.Load(0);

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function Clear(T)

if T.read, error('read only'); end

T.Clear2(0);

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function ClearFields(T)

if T.read, error('read only'); end

if T.State == 0
    T.Clear2(0);
else
    T.Clear2(1);
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function ClearData(T)

if T.read, error('read only'); end
if T.State < 2, error('no fields'); end

T.Clear2(2);

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = private)
function Clear2(T, state)

if state < 1
    T.digs      = [];
    T.blockSize = [];
end

if state < 2
    T.fields  = {};
    T.classes = {};
    T.sizes   = {};
    T.auto    = [];
end

T.parts = 0;
T.c     = 0;
T.mod   = false;
T.cnt   = inf;
T.buf   = [];
T.n     = 0;

T.SetBlockSize;

T.DeleteBlocks(state < 2);

if state > 0
    T.SaveHeader;
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Hidden)
function U = saveobj(T)

T.Load(0);

U = T;

U.opts = intersect(U.opts, 'r');

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function delete(T)

T.Load(0);

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function Copy(T, b, over)

if nargin < 3, over = false; end

if T.read, error('read only'); end
if (T.State > 0) && ~over, error('array is already initialized, overwrite not specified'); end

if b.State == 0
    T.Clear;
else
    T.Init(b.BlockSize, b.Digits);
    if b.State >= 2
        T.Append(b);
        T.Flush;
    end
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function state = State(T)

if isempty(T.digs)
    state = 0;
elseif isempty(T.fields)
    state = 1;
elseif T.n == 0
    state = 2;
else
    state = 3;
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function digs = Digits(T)

digs = T.digs;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function [blockSize, blockSize2] = BlockSize(T)

blockSize  = T.blockSize;
blockSize2 = T.blockSize2;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function nf = NumFields(T)

nf = numel(T.fields);

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function f = Fields(T, i)

if nargin < 2
    f = T.fields;
else
    f = T.fields{i};
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function c = Classes(T, i)

if nargin < 2
    c = T.classes;
elseif ischar(i)
    [ans, j] = ismember(i, T.fields);
    if j == 0, error('invalid field "%s"', i); end
    c = T.classes{j};
else
    c = T.classes{i};
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function s = Sizes(T, i)

if nargin < 2
    s = T.sizes;
elseif ischar(i)
    [ans, j] = ismember(i, T.fields);
    if j == 0, error('invalid field "%s"', i); end
    s = T.sizes{j};
else
    s = T.sizes{i};
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = protected)
function Display(T, path)

fprintf('\n');

if iscellstr(path)
    path = sprintf('{''%s'', ''%s''}', path{1}, path{2});
else
	path = sprintf('''%s''', path);
end

if isempty(T.blockSize)
    bs = '[]';
elseif T.blockSize > 0
    bs = sprintf('%u items', T.blockSize);
elseif isempty(T.blockSize2)
    bs = sprintf('%u bytes', -T.blockSize);
else
    bs = sprintf('%u bytes (%u items)', -T.blockSize, T.blockSize2);
end

fprintf('     CLASS: ''%s''\n', class(T));
fprintf('      PATH: %s\n', path);
fprintf('BLOCK SIZE: %s\n', bs);
fprintf('    BLOCKS: %u\n', T.parts);
fprintf('     ITEMS: %u\n', T.n);
fprintf('\n');

fprintf('FIELDS:\n');
if T.auto
    names = {['*' T.fields{1}]};
else
    names = T.fields;
end
nlen = max(cellfun(@numel, names));
for i = 1 : numel(T.fields)
    siz = sprintf('%ux', T.sizes{i});
    fprintf('  %*s: %sN %s\n', nlen, names{i}, siz, T.classes{i});
end
fprintf('\n');

fprintf('ACCESS FLAGS:\n');
s.readOnly = sprintf('%u', T.read);
s.verbose  = sprintf('%u', T.verb);
if isempty(T.perm)
    s.transpose = '[]';
elseif islogical(T.perm)
    s.transpose = sprintf('%u', T.perm);
else
    s.transpose = mat2str(T.perm);
end
names = fieldnames(s);
nlen = max(cellfun(@numel, names));
for i = 1 : numel(names)
    name = names{i};
    fprintf('  %*s: %s\n', nlen, name, s.(name));
end
fprintf('\n');

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function ni = NumItems(T)

ni = T.n;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function d = Read(T, inds, fields)

if T.State < 2, error('no fields'); end

if islogical(inds)
    inds = find(inds);
    back = [];
elseif isempty(inds)
    back = [];
elseif any(diff(inds) < 1)
    [inds, ans, back] = unique(inds);
else
    back = [];
end
if ~isempty(inds)
    if (inds(1) < 1) || (inds(end) > T.n), error('invalid indices'); end
end

if nargin < 3
    fields = T.fields;
    auto   = T.auto;
    is     = 1 : numel(T.fields);
else
    if ischar(fields)
        fields = {fields};
        auto = true;
    else
        if isempty(fields), error('field mismatch'); end
        auto = false;
    end
    [ans, is] = ismember(fields, T.fields);
    if any(is == 0), error('field mismatch'); end
end

d = struct;
for i = is
    d.(T.fields{i}) = repmat(Filler(T.classes{i}), [T.sizes{i}, numel(inds)]);
end

if ~isempty(inds)
    for c = 1 : T.parts
        match = (floor((inds - 1) / T.blockSize2) + 1 == c);
        pinds = inds(match) - (c - 1) * T.blockSize2;
        if isempty(pinds), continue; end
        off = find(match, 1, 'first') - 1;
        m = T.Buffer(c);
        for i = is
            colons = T.Colons(i);
            d.(T.fields{i})(colons{:}, off + 1 : off + numel(pinds)) = m.Data.(T.fields{i})(colons{:}, pinds);
        end
    end
end

if ~isempty(back)
    for i = is
        colons = T.Colons(i);
        d.(T.fields{i}) = d.(T.fields{i})(colons{:}, back);
    end
end

if auto, d = d.(fields{1}); end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = private)
function m = Buffer(T, c)

if c == T.c
    m.Data = T.buf;
else
    m = T.BlockMap(c);
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = protected)
function m = BlockMap(T, c)

T.Load(c);

m.Data = T.buf;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function nb = NumBlocks(T)

nb = T.parts;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function d = ReadBlock(T, c, varargin)

if T.State < 2, error('no fields'); end

if ~isscalar(c), error('invalid block number'); end
if (c < 1) || (c > T.parts), error('invalid block number'); end

ind1 = (c - 1) * T.blockSize2 + 1;
ind2 = min(c * T.blockSize2, T.n);

d = T.Read(ind1 : ind2, varargin{:});

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function prev = Transpose(T, perm)

prev = T.perm;

if nargin < 2, return; end

% TODO: check
T.perm = perm;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Hidden)
function varargout = subsref(T, ids)

switch ids(1).type
case '.'

    [varargout{1 : nargout}] = builtin('subsref', T, ids);

case '()'

    if T.State < 2, error('no fields'); end
    if ~T.auto, error('direct indexing not available for arrays that use field names'); end
    nd = numel(T.sizes{1}) + 1;

    if isempty(T.perm)
        perm = 1 : nd;
    elseif islogical(T.perm)
        if T.perm
            perm = [2 1];
        else
            perm = [1 2];
        end
    else
        perm = T.perm;
    end
    if numel(perm) ~= nd, error('wrong number of transpose dimensions'); end
    [ans, last] = max(perm);

    if numel(ids) > 1, error('invalid reference'); end
    ni = numel(ids.subs);
    if (ni < last) || (ni > nd), error('wrong number of dimensions'); end
    if (ni == last) && (ni < nd), error('wrong number of dimensions'); end

    for i = 1 : ni
        if i == last
            if isequal(ids.subs{i}, ':'), error('must supply indices for dimension %u', last); end
        else
            if ~isequal(ids.subs{i}, ':'), error('can only supply indices for dimension %u', last); end
        end
    end

    d = permute(T.Read(ids.subs{last}), perm);

    if ni < nd
        sizes = size(d);
        sizes(end + 1 : nd) = 1;
        sizes = [sizes(1 : ni - 1), prod(sizes(ni : end))];
        d = reshape(d, sizes);
    end

    varargout{1} = d;

otherwise

    error('invalid reference');

end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function BlockOp(T, f, write)

if nargin < 3, write = false; end

if write && T.read, error('read only'); end
if T.State < 2, error('no fields'); end

for c = 1 : T.parts

    T.Load(c);

    buf = struct;
    for i = 1 : numel(T.fields)
        colons = T.Colons(i);
        buf.(T.fields{i}) = T.buf.(T.fields{i})(colons{:}, 1 : T.cnt);
    end
    if T.auto, buf = buf.v; end

    if write

        buf = f(buf);

        [buf, cnt] = T.CheckData(buf);
        if cnt ~= T.cnt, error('cannot change the number of items'); end

        T.buf = buf;
        for i = 1 : numel(T.fields)
            colons = T.Colons(i);
            T.buf.(T.fields{i})(colons{:}, T.cnt + 1 : T.blockSize2) = Filler(T.classes{i});
        end

        T.mod = true;

    else

        f(buf);

    end

end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Static)
function a = InnerProd(T, U, f1, f2)

if nargin < 3, f1 = 'v'; end
if nargin < 4, f2 = f1 ; end

cls = T.Classes(f1);
if ~strcmp(cls, U.Classes(f2)), error('class mismatch'); end
if ~isscalar(T.Sizes(f1)) || ~isscalar(U.Sizes(f2)), error('matrices required'); end

a = zeros(T.NumItems, U.NumItems, cls);

for i1 = 1 : T.BlockSize : T.NumItems

    i2 = min(i1 + T.BlockSize - 1, T.NumItems);
    t = T.Read(i1 : i2, f1)';

    for j1 = 1 : U.BlockSize : U.NumItems

        j2 = min(j1 + U.BlockSize - 1, U.NumItems);
        u = U.Read(j1 : j2, f2);

        a(i1 : i2, j1 : j2) = t * u;

    end
    
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = private)
function SaveHeader(T)

s.type      = class(T);
s.blockSize = T.blockSize;
s.fields    = T.fields;
s.classes   = T.classes;
s.sizes     = T.sizes;
s.auto      = T.auto;

T.SaveBlock(0, s);

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = private)
function Load(T, c)

if c == T.c, return; end

if T.mod
    s = struct;
    s.num = T.cnt;
    s.var = struct;
    for i = 1 : numel(T.fields)
        colons = T.Colons(i);
        s.var.(T.fields{i}) = T.buf.(T.fields{i})(colons{:}, 1 : T.cnt);
    end
    T.SaveBlock(T.c, s);
end

T.c   = c;
T.mod = false;

if T.c == 0
    T.cnt = inf;
    T.buf = [];
elseif T.c <= T.parts
    s = T.LoadBlock(T.c);
    T.cnt = s.num;
    T.buf = struct;
    for i = 1 : numel(T.fields)
        colons = T.Colons(i);
        T.buf.(T.fields{i}) = s.var.(T.fields{i});
        if ~T.read
            T.buf.(T.fields{i})(colons{:}, T.cnt + 1 : T.blockSize2) = Filler(T.classes{i});
        end
    end
else
    if numel(num2str(T.c)) > T.digs, error('maximum number of parts exceeded'); end
    T.cnt = 0;
    T.buf = struct;
    for i = 1 : numel(T.fields)
        T.buf.(T.fields{i}) = repmat(Filler(T.classes{i}), [T.sizes{i}, T.blockSize2]);
    end
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = private)
function c = Colons(T, i)

c = repmat({':'}, 1, numel(T.sizes{i}));

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Abstract, Access = protected)
	[list, digs] = FindBlocks(T);
    DeleteBlocks(T, head);
	s = LoadBlock(T, c, field);
    SaveBlock(T, c, s);
end

%-----------------------------------------------------------------------------------------------------------------------

end

%***********************************************************************************************************************

function a = Null(a)

end

%***********************************************************************************************************************

function a = Filler(cls)

if strcmp(cls, 'cell')
    a = cell(1);
else
    a = zeros(1, cls);
end

end

%***********************************************************************************************************************

function b = Bytes(cls)

switch cls
case 'char'
    a = 'x'; %#ok
case 'logical'
    a = true; %#ok
otherwise
    try
        a = zeros(1, cls); %#ok
    catch
        a = []; %#ok
    end
end

w = whos('a');
b = w.bytes;

end