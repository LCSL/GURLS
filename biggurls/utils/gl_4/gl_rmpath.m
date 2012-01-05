function gl_rmpath(path)

if nargin < 1
    RemovePath(cd);
else
    current = cd(path);
    RemovePath(cd);
    cd(current);
end

return;

%***********************************************************************************************************************

function RemovePath(path)

files = dir(path);

for i = numel(files) : -1 : 1
    if ~files(i).isdir                                            , continue; end
    if files(i).name(1) == '.'                                    , continue; end
    if ismember(files(i).name, {'private'})                       , continue; end
    if exist(fullfile(path, files(i).name, '.nopath'    ), 'file'), continue; end
    if exist(fullfile(path, files(i).name, '.nopath.txt'), 'file'), continue; end
    RemovePath(fullfile(path, files(i).name));
end

warning('off', 'MATLAB:rmpath:DirNotFound');
rmpath(path);
warning('on', 'MATLAB:rmpath:DirNotFound');

return;