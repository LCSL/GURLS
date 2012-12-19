classdef bigarray_mem < bigarray

%-----------------------------------------------------------------------------------------------------------------------

properties (Access = private)
    store;
    storeDigs;
end

%-----------------------------------------------------------------------------------------------------------------------

methods
function T = bigarray_mem(varargin)

T@bigarray(varargin{:});

T.store     = {};
T.storeDigs = [];

T.Construct;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Static, Hidden)
function T = loadobj(U)

T = U;

T.Construct;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Hidden)
function display(T)

T.Display('');

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = protected)
function [list, digs] = FindBlocks(T)

list = 0 : numel(T.store) - 1;
digs = T.storeDigs;

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = protected)
function DeleteBlocks(T, head)

if head
    T.store     = {};
    T.storeDigs = [];
else
    T.store(2 : end) = [];
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = protected)
function s = LoadBlock(T, c, field)

if nargin < 3, field = ''; end

s = T.store{c + 1};

if ~isempty(field), s = s.(field); end

end
end

%-----------------------------------------------------------------------------------------------------------------------

methods (Access = protected)
function SaveBlock(T, c, s)

T.store{c + 1} = s;

if c == 0
    T.storeDigs = T.digs;
end

end
end

%-----------------------------------------------------------------------------------------------------------------------

end