function path = biggurls_install()

		basePath = biggurls_root;

		run(fullfile(basePath,'../gurls/utils/gurls_install'	));
		run(fullfile(basePath,'utils/adm_install'				));

		addpath(basePath										 );
		addpath(fullfile(basePath,'kernels'						));
		addpath(fullfile(basePath,'optimizers'					));
		addpath(fullfile(basePath,'paramsel'					));
		addpath(fullfile(basePath,'pred'						));
		addpath(fullfile(basePath,'perf'						));
		addpath(fullfile(basePath,'norm'						));
		addpath(fullfile(basePath,'utils'						));
		addpath(fullfile(basePath,'summary'						));
		addpath(fullfile(basePath,'confidence'					));
		addpath(fullfile(basePath,'demo'						));
		addpath(fullfile(basePath,'distributed'					));
		addpath(fullfile(basePath,'demo/code'					));
		addpath(fullfile(basePath,'adm/bin'						));
		addpath(fullfile(basePath,'adm/utils'					));
		addpath(fullfile(basePath,'adm/'						));
		addpath(fullfile(basePath,'utils/bigarray_1'			));
		addpath(fullfile(basePath,'utils/gl_4'					));
