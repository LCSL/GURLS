% test ISTA on elastic net
name='ExampleExperiment';
opt = gurls_defopt(name);
opt.seq={'split:ho','paramsel:hoista','rls:ista','pred:primal','perf:macroavg'};
opt.process{1}=[2,2,2,0,0]; 
opt.process{2}=[3,3,3,2,2]; 
%opt.INSTASparsitylvl=1;
opt.INSTASAccuReq = 0.3;
gurls(Xtr,ytr,opt,1);
gurls(Xte,yte,opt,2);