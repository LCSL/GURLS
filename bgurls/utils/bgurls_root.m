function path = bgurls_root

[a,b,c] = fileparts(mfilename('fullpath'));
[a,b,c] = fileparts(a);
path = a;
