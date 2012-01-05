function store = gl_savepath(use, store, varargin)

if ~iscell(store), return; end

test = gl_abspath(store, varargin{:});

if ~strcmp(test, use), store = use; end
% TODO: attempt to fixup relative path instead?

return;