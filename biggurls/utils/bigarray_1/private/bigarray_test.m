function bigarray_test

d    = 10000;
n    = 105;
bsiz = 13;
base = fullfile(workdir, 'bmtest');

% Generate fake data blocks.
for i = 1 : ceil(n / bsiz)
    a1{i} = rand(d, min(n - bsiz * (i - 1), bsiz), 'double');
end

% Compute inner product using ordinary MATLAB.
b1 = cat(2, a1{:});
b1 = b1' * b1;

% Now dump the data blocks into a bigarray.
a2 = bigarray_bin(base, 'v');
a2.Init(bsiz);
for i = 1 : numel(a1)
    a2.Append(a1{i});
end

% Compute inner product using the bigarray (then delete the array).
b2 = bigarray.InnerProd(a2, a2);
a2.Clear;

% See if the results are the same.
hist(abs(b2(:) - b1(:)), 100);

return;