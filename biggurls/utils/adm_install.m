function path = adm_install()

		basePath = biggurls_root;
		basePath = fullfile(basePath,'adm');
		
		srcdir 	= fullfile(basePath, 'src');
		bindir 	= fullfile(basePath, 'bin');
		
		mex('-outdir', bindir, fullfile(srcdir, 'admSetup.cpp'));
		mex('-outdir', bindir, fullfile(srcdir, 'admDismiss.cpp'));
		mex('-outdir', bindir, 	fullfile(srcdir, 'getWork.cpp'), ...
								fullfile(srcdir, 'lockManagement.cpp'),... 
								fullfile(srcdir, 'workManagement.cpp'));
		mex('-outdir', bindir, 	fullfile(srcdir, 'reportWork.cpp'), ...
								fullfile(srcdir, 'lockManagement.cpp'),... 
								fullfile(srcdir, 'workManagement.cpp'));
		
		addpath(bindir);
