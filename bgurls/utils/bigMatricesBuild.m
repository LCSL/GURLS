function [] = bigMatricesBuild(files)
% Pre-computes and stores matrices XtX, Xty, XvatXva, Xvatyva
% INPUTS:
% -files = structure containing the fields:
%     -Xtrain_filename = prefix of files that make the bigarray for input
%       training data
%     -ytrain_filename = prefix of files that make the bigarray for ouput
%       training data
%     -Xva_filename = prefix of files that make the bigarray for input
%       validation data
%     -ytva_filename = prefix of files that make the bigarray for ouput
%       validation data
%     -XtX_filename = name of file where matric X'*X will be stored
%     -Xty_filename = name of file where matric X'*y will be stored
%     -XvatXva_filename = name of file where matric Xva'*Xva will be stored
%     -Xvatyva_filename = name of file where matric Xva'*yva will be stored

%% Set up distributed matrix-matrix multiplications (XtX, Xty, XvatXva,Xvatyva) with gdm

X = bigarray.Obj(files.Xtrain_filename);
y = bigarray.Obj(files.ytrain_filename);

admMatMultPrepare(X,X,files.XtX_filename);
admMatMultPrepare(X,y,files.Xty_filename); 

Xva = bigarray.Obj(files.Xva_filename);
yva = bigarray.Obj(files.yva_filename);

admMatMultPrepare(Xva,Xva,files.XvatXva_filename);
admMatMultPrepare(Xva,yva,files.Xvatyva_filename);

%% Actually run distributed matrix-matrix multiplications with gdm (run this on multiple machines after you see the message!!).
% builds matrices XtX, Xty, XvatXva, Xvatyva

fprintf('Run this\n\tbigMatricesBuild(''<FILES>'')\non as many machines as you please\n');
admMatMultRun(files.XtX_filename);
admMatMultRun(files.Xty_filename); 
admMatMultRun(files.XvatXva_filename);
admMatMultRun(files.Xvatyva_filename);
