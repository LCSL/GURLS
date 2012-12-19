function path = bgurls_install()

        basePath = bgurls_root;

        run(fullfile(basePath,'../gurls/utils/gurls_install'	));
        run(fullfile(basePath,'utils/adm_install'				));

        addpath(basePath										 );
        addpath(fullfile(basePath,'optimizers'					));
        addpath(fullfile(basePath,'paramsel'					));
        addpath(fullfile(basePath,'pred'						));
        addpath(fullfile(basePath,'perf'						));
        addpath(fullfile(basePath,'utils'						));
        addpath(fullfile(basePath,'demo'						));
        addpath(fullfile(basePath,'distributed'					));
        addpath(fullfile(basePath,'demo/code'					));
        addpath(fullfile(basePath,'adm/bin'						));
        addpath(fullfile(basePath,'adm/utils'					));
        addpath(fullfile(basePath,'adm/'						));
        addpath(fullfile(basePath,'utils/bigarray_1'			));
        addpath(fullfile(basePath,'utils/gl_4'					));
