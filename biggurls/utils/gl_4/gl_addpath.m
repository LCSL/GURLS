function gl_addpath(path)

if nargin < 1
    SetPath(cd);
else
    current = cd(path);
    SetPath(cd);
    cd(current);
end

return;

%***********************************************************************************************************************

function SetPath(path)

files = dir(path);

for i = numel(files) : -1 : 1
    if ~files(i).isdir                                            , continue; end
    if files(i).name(1) == '.'                                    , continue; end
    if ismember(files(i).name, {'private'})                       , continue; end
    if exist(fullfile(path, files(i).name, '.nopath'    ), 'file'), continue; end
    if exist(fullfile(path, files(i).name, '.nopath.txt'), 'file'), continue; end
    SetPath(fullfile(path, files(i).name));
end

addpath(path);

return;