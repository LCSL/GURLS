function path = gurls_install()

basePath = gurls_root();
addpath(basePath							);
addpath(fullfile(basePath,'kernels'		)	);
addpath(fullfile(basePath,'optimizers'	)	);
addpath(fullfile(basePath,'paramsel'	)	);
addpath(fullfile(basePath,'pred'		)	);
addpath(fullfile(basePath,'perf'		)	);
addpath(fullfile(basePath,'norm'		)	);
addpath(fullfile(basePath,'utils/'		)	);
addpath(fullfile(basePath,'summary'		)	);
addpath(fullfile(basePath,'confidence'	)	);

