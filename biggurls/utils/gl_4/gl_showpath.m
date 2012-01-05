function gl_showpath

list = strread(path, '%s', 'delimiter', pathsep);
root = matlabroot;

fprintf('\n');

for i = 1 : numel(list)
    if isempty(strmatch(root, list{i}))
        fprintf('        %s\n', list{i});
    end
end

fprintf('\n');

return;