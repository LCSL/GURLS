name = 'ExampleExperiment'
opt = defopt(name);
opt.seq = ...
{'split:ho','paramsel:bfprimal','rls:landweberprimal_lgs', ...
 'pred:primal','perf:macroavg'};
%{'split:ho','kernel:rbf','paramsel:bfdual','rls:landweberdual', ...
% 'predkernel:traintest','pred:dual','perf:macroavg'};
opt.process{1} = [2,2,2,0,0];
opt.process{2} = [3,3,3,2,2];
opt.nholdouts = 5;
opt.paramsel.guesses = linspace( 100, 1000, 301 );
opt.paramsel.optimizer = @rls_landweberprimal_lgs;
opt.paramsel.sigma = 1;
gurls (Xtrain, Ytrain, opt, 1)
gurls (Xtest, Ytest, opt, 2)