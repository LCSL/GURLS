% installing 'libsvmread' from LIBSVM
if (exist('libsvmread') == 0)
    cd('/Users/Jeremiah/GitHub/GURLS/gurls/optimizers/mkl/data/libsvm');
    mex CFLAGS="\$CFLAGS -std=c99" -largeArrayDims libsvmread.c
end

% find _scale.txt files

file_list = dir('*_scale.txt');
file_name = strrep({file_list.name}, '_scale.txt', '');

% convert
for idx = 1:length(file_name)
    [y, X] = libsvmread(['./', file_name{idx}, '_scale.txt']);
    csvwrite(['../', file_name{idx}, '.csv'], full([y, X]));
end