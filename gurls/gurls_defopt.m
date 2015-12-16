function opt = gurls_defopt(expname)

    if (nargin < 1)
        fprintf('Your experiment must be given a unique name\n');
    end

    opt = GurlsOptions();
    
    opt.newprop('time', {});

    %% Experiment options
    opt.newprop('name', expname); % We can make this argument mandatory; Unique name per experiment
    opt.newprop('savefile', [opt.name '.mat']);
    opt.newprop('save', true);

    %% Algorithm options
    opt.newprop( 'kernel', struct());
    opt.newprop( 'predkernel', struct());
    opt.kernel.type = 'rbf';
    opt.newprop( 'singlelambda', @median); % give the function for combining your lambdas; could be max, min or mean, for instance.
    opt.newprop( 'smallnumber' ,1e-8); % lambda is searched between [min(eig_r, opt.smallnumber), eig_1], where r = rank, eig_1 = max eig val.

    %% Random SVD options
    opt.newprop('eig_percentage', 5); %percentage of eigenvectors to be used in the randomized SVD

    %% Iteartive RLS options
    opt.newprop('IterRLSMaxIter', 1000);
    opt.newprop('IterRLSMinIter', 5);
    opt.newprop('IterRLSStopTol', 0.001);
    opt.newprop('IterRLSSeriesType','geometric');
    
    %% GD Options
    opt.newprop('gd', struct());
    opt.gd.method = 0; % standard gradient descent
    opt.gd.maxiter = 1000;
    opt.gd.singleiter = @median;
    opt.gd.eta_numerator = 1;
    opt.gd.nu = 1;

    %% CG Options
    opt.newprop('cg', struct());
    opt.cg.maxiter = 1000;
    opt.cg.singleiter = @median;

    %% Output options
    opt.newprop('hoperf', @perf_macroavg);
    opt.newprop('nholdouts', 1);
    
    %% Data option
    opt.newprop( 'hoproportion', 0.2);
    opt.newprop( 'nlambda', 20);
    opt.newprop( 'nsigma', 25);
    
    %% Quiet
    % Currenty either 0 or 1; levels of verbosity may be implemented later;
    opt.newprop( 'verbose', 1);

    %% Version info
    opt.newprop( 'version', 2.0);

    %% Online
    opt.newprop( 'epochs', 4);
    opt.newprop( 'subsize', 50);
    opt.newprop( 'calibfile', 'foo');

    %% Random features options
    opt.newprop( 'randfeats', struct());
    opt.randfeats.D = 500;
    opt.randfeats.samplesize = 100;

    %% Nystrom options
    opt.newprop( 'nystrom', struct());
    opt.nystrom.m = 100;
    
    %% Preprocessing / Dimensionality Reduction options
    opt.newprop( 'preproc', struct() );
    opt.preproc.kernel = struct;
    opt.preproc.kernel.kernel = 'linear';
    opt.preproc.center_data = true;
    opt.preproc.n_dims = 5;
    
    %% opt base
    opt.newprop( 'jobid', 1);
    opt.newprop( 'seq', {});
    opt.newprop( 'process', {});
end
