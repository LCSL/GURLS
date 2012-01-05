function [] = bgTrainRun(wpath)
		admMatMultRun(fullfile(wpath,'XtX.mat'));
		admMatMultRun(fullfile(wpath,'Xty.mat')); 
		admMatMultRun(fullfile(wpath,'XvatXva.mat'));
		admMatMultRun(fullfile(wpath,'Xvatyva.mat'));
