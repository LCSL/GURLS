load('./data/yeast_data.mat');

for i = 1:5
	name = ['linloo_' num2str(i)];
	opt = defopt(name);
	opt.seq = {'paramsel:loocvprimal','rls:primal','pred:primal','perf:precrec','perf:macroavg'};
	opt.process{1} = [2,2,0,0,0];
	opt.process{2} = [3,3,2,2,2];
	gurls (Xtr, ytr, opt,1)
	gurls (Xte, yte, opt,2)
	
	name = ['rbfloo_' num2str(i)];
	opt = defopt(name);
	opt.seq = {'paramsel:siglam', 'kernel:rbf', 'rls:dual', 'predkernel:traintest', 'pred:dual', 'perf:macroavg', 'perf:precrec'};
	opt.process{1} = [2,2,2,0,0,0,0];
	opt.process{2} = [3,3,3,2,2,2,2];
	gurls (Xtr, ytr, opt,1)
	gurls (Xte, yte, opt,2)
	
	
	name = ['linho_' num2str(i)];
	opt = defopt(name);
	opt.seq = {'split:ho','paramsel:hoprimal','rls:primal','pred:primal','perf:macroavg','perf:precrec'};
	opt.process{1} = [2,2,2,0,0,0];
	opt.process{2} = [3,3,3,2,2,2];
	gurls (Xtr, ytr, opt,1)
	gurls (Xte, yte, opt,2)
	
	
	name = ['rbfho_' num2str(i)];
	opt = defopt(name);
	opt.seq = {'split:ho', 'paramsel:siglamho', 'kernel:rbf', 'rls:dual', 'predkernel:traintest', 'pred:dual', 'perf:macroavg', 'perf:precrec'};
	opt.process{1} = [2,2,2,2,0,0,0,0];
	opt.process{2} = [3,3,3,3,2,2,2,2];
	gurls (Xtr, ytr, opt,1)
	gurls (Xte, yte, opt,2)

end

filestr = {'linloo','rbfloo','linho','rbfho'};
nRuns = {5,5,5,5};
fields = {'perf.ap','perf.acc'};
plotopt.titles = {'Model Comparison - Accuracy','Model Comparison - Precision'};
plotopt.class_names = {'CYT','NUC','MIT','ME3'};
summary_plot(filestr,fields,nRuns,plotopt)  
summary_overall_plot(filestr,fields,nRuns,plotopt)  
plot_times(filestr,nRuns)

