function model = train(X, y, varargin)
    opt = prodOptions(varargin);
    y = analyzeData(X, y, opt);
    chooseMatches(opt);
    chooseCVAlg(opt); 
    modelselection(opt);
    
    seq = {};
    if isequal(opt.partuning,'ho') && ~isequal(opt.pars, 'none')
        seq = {'split:ho'};
    end
    
    switch opt.pars
        case 'none'
            if ~isequal(opt.algorithm,'lrls') 
               seq = [seq, {['kernel:' opt.kernelfun]}];
            end
        case 'reg'
            if isequal(opt.algorithm,'lrls')
                seq = [seq, {['paramsel:' opt.cvalgo]}];
            else
               seq = [seq, {['kernel:' opt.kernelfun], ['paramsel:' opt.cvalgo]}];
            end
        case {'ker', 'all'}
            seq = [seq, {['paramsel:' opt.cvalgo], ['kernel:' opt.kernelfun]}];
    end
    
    model = GurlsOptions(opt);
    
	opt.seq = [seq, {['rls:' opt.algname]}];
	opt.process{1} = 2*ones(1, numel(opt.seq));
	gurls(X, y, opt, 1);
    
    model.newprop('paramsel', opt.paramsel);
    model.newprop('rls', opt.rls);
    model.newprop('kernel', opt.setting.kernel);
    
    model.seq = {};
    model.process{1} = [];
    if ~(isequal(opt.algorithm,'lrls') || (isequal(opt.kernel.type, 'linear') && isequal(opt.algname, 'dual')))
         model.seq = {'predkernel:traintest'};
         model.process{1} = 1;
    end
    model.seq = [model.seq,  {['pred:' opt.predname], ['perf:' opt.perfm]}];
    model.process{1} = [model.process{1}, 1, 1];
end

function opt = prodOptions(lst)
fprintf('\n\n\nGurls 2.0\n');

    if numel(lst) == 1
        opt = lst{1};
        if ~isa(opt, 'GurlsOptions')
            opt = GurlsOptions(opt);
        end
    elseif ~mod(numel(lst),2)
        opt = defopt('');
        for i=1:numel(lst)/2
            if ~ischar(lst{2*i-1})
                error('The %d-th option name is not a string', i);
            elseif ~isequal(regexp(lst{2*i-1}, '\w*(\.\w*)*'),1)
                error('The option "%s" is in the wrong format', lst{2*i-1});
            else
                opt.newprop(lst{2*i-1}, lst{2*i});
            end
        end
    else
        error('The option list is wrongly formatted');
    end
end


function selectDatatype(X, opt)
    ok = false;
    if isprop(opt,'datatype')
        switch opt.datatype
            case 'vector'
                if ~isnumeric(X)
                    disp('The supplied input data is not a matrix');
                elseif issymmetric(X)
                    disp('The supplied input matrix is a matrix, but it is symmetric.');
                else
                    ok = true;
                end
            case 'kernel'
                if ~isnumeric(X)
                    disp('The supplied input data is not a matrix, thus it cannot be of type "kernel"');
                elseif ~issymmetric(X)
                    disp('The supplied input matrix is not symmetric, thus the chosen datatype cannot be "kernel"');
                else
                    ok = true;
                end
            case 'string'
                error('The option "datatype" can be "vector" or "kernel"');
%                 if ~iscellstr(X)
%                     disp('The supplied input data is not of a cell of strings, thus it cannot be of type "string".');
%                 else
%                     ok = true;
%                 end
            case 'user-defined'
                error('The option "datatype" can be "vector" or "kernel"');
                %ok = true;
            otherwise
                error('The option "datatype" can be "vector" or "kernel"');
                %opt.datatype = unknelem('datatype', opt.datatype, '(vector|kernel|string|user-defined)', '');
                %ok = false;
        end
    end
    
    if ~isprop(opt,'datatype') || ~ok
        opt.newprop('datatype', []);
        if isnumeric(X)
                opt.datatype = 'vector';
                [a, b] = size(X);
                if (a == b) && issymmetric(X)
                    opt.datatype = 'kernel';
                end

%         elseif iscellstr(X)
%             opt.datatype = 'string';
%         elseif iscell(X)
%             opt.datatype = 'user-defined';
        else
            error('Input data not recognizable');
        end
    end
    if ~isprop(opt,'datatype')
        fprintf('No datatype has been supplied by the user. ');
    end
    fprintf('The selected datatype for the input is "%s"\n', opt.datatype);
end

function y = analyzeData(X, y, opt)
    
    selectDatatype(X, opt);
                
    opt.newprop('n', -1);
    opt.newprop('d', -1);
    switch opt.datatype
        case 'vector'
            [opt.n, opt.d] = size(X);
        case 'kernel'
            opt.n = size(X,1);
        case 'string'
            opt.n = numel(X);
    end
    
    
    [a, T] = size(y);
    if a ~= opt.n && T == opt.n
        y = y';
        T = a;
    elseif a ~= opt.n && T ~= opt.n
        error('There are %d examples but %d labels', opt.n, a);
    end
            
    
    if ~isprop(opt,'problem') || ~ischar(opt.problem) || ...
        ~isequal(regexp(opt.problem, '(regression|classification)'), 1)
        opt.newprop('problem', 'classification');
        if max(max(abs(double(int32(y))-y))) > eps
            opt.problem ='regression'; 
        end
        fprintf('The problem has been set to %s.\n', opt.problem);
    end
    
    if isequal(opt.problem, 'classification') && max(max(abs(double(int32(y))-y))) < eps
        if T==1
            ind = sort(unique(floor(y + 0.5)));
            opt.newprop('labeldict', ind);
            if numel(ind) > 2
                opt.newprop('yformat', 'compact-multiclass');
                yp = zeros(numel(y),1);
                for i=1:size(ind)
                    yp(y==ind(i)) = i;
                end
                I = 2*eye(size(ind,1)) - 1;
                y = I(yp,:); 
                T = size(ind);
            else
                opt.newprop('yformat', 'compact-twoclass');
                z = -ones(size(y,1),1);
                z(y==max(y)) = 1;
                y = z;
                T = 2;
            end
        else
            opt.newprop('yformat', 'extended');
        end
        opt.newprop('numoutputs', T);
    elseif isequal(opt.problem, 'classification') && max(max(abs(double(int32(y))-y))) > eps
       warning(['You selected a classification problem but the labels seem to be reals.'...
                'Have you considered to perform a regression problem instead?']);
       opt.newprop('yformat', 'notinteger');
    elseif isequal(opt.problem, 'regression')
        if max(max(abs(double(int32(y))-y))) < eps
            warning(['You selected a regression problem but the labels seem to be integers.'...
                'Have you considered to perform a classification problem instead?']);
        end
        opt.newprop('numoutputs', size(y,2));
    end
    
    opt.newprop('setting', struct());
    opt.setting.datatype = opt.datatype;
    opt.setting.n = opt.n;
    opt.setting.d = opt.d;
    opt.setting.problem = opt.problem;
    opt.setting.numoutputs = opt.numoutputs;
    if isprop(opt, 'labeldict')
        opt.setting.labeldict = opt.labeldict;
    end
    if isprop(opt, 'yformat')
        opt.setting.yformat = opt.yformat;
    end
    
end


function kernel = selectKernel(opt, isLRLS)
    ok = false;
    
    if ~isprop(opt, 'kernelfun') && ischar(opt.kernel)
        opt.newprop('kernelfun', opt.kernel);
        opt.kernel = struct();
    end
    
    if isprop(opt, 'kernelfun')
        if ischar(opt.kernelfun)
            if ~isequal('', unknelem('kernel', opt.kernelfun, '(linear|datatype|user-supplied|rbf|chisquared|quasiperiodic|string)', ''))
                switch opt.datatype
                    case 'vector'
                        if isLRLS
                            kernel = allowelem('kernel', 'vector datatype and lrls algorithm', opt.kernelfun, 'linear', 'linear');
                        else
                            kernel = allowelem('kernel', 'vector datatype', opt.kernelfun, '(linear|rbf|chisquared|quasiperiodic)', 'rbf');
                        end
                    case 'kernel'
                        kernel = allowelem('kernel', 'kernel datatype', opt.kernelfun, 'datatype', 'datatype');
                    case 'string'
                        kernel = allowelem('kernel', 'string datatype', opt.kernelfun, 'string', 'string');
                    case 'user-defined'
                        error('Kernel for the user-defined datatype must be explicitly provided');
                end
                ok = true;
            end
        elseif isa(opt.kernelfun, 'function_handle')
            kernel = 'user-supplied';
        else
            error('Unrecognizable option "kernelfun"');
        end
    end
    
    if ~isprop(opt, 'kernelfun') || ~ok
        switch opt.datatype
            case 'vector'
                if isLRLS
                    kernel = 'linear';
                else
                    kernel = 'rbf';
                end
            case 'kernel'
                kernel = 'datatype';
            case 'string'
                kernel = 'string';
            case 'user-defined'
                error('Kernel for unknown datatype must be explicitly provided');
        end
        if ~ok
            fprintf('The selected kernel is "%s".\n', kernel);
        end
        opt.newprop('kernelfun', kernel);
    end
end

function algorithm = selectAlgo(opt)
    ok = false;
    algorithm = 'user-supplied';
    
    if isprop(opt, 'algorithm')
        if ischar(opt.algorithm)
            if ~isequal('', unknelem('algorithm', opt.algorithm, '(lrls|krls|krlsrf|gpr)', ''))
                switch opt.datatype
                    case 'vector'
                        algorithm = allowelem('algorithm', 'vector datatype', opt.algorithm, '(lrls|krls|krlsrf|gpr)', 'krls');
                    case 'kernel'
                        algorithm = allowelem('algorithm', 'kernel datatype', opt.algorithm, '(krls|krlsrf|gpr)', 'krls');
                    case {'string', 'user-defined'}
                        algorithm = allowelem('algorithm', 'non-vector datatypes', opt.algorithm, '(krls|krlsrf|gpr)', 'krls');
                end
                ok = true;
            end
        elseif isa(opt.algorithm, 'function_handle')
%             algorithm = 'user-supplied';
%             opt.newprop('algofun', opt.algorithm);
%             opt.setting.algofun = opt.algofun;
%             ok = true;
            error('The option "algorithm" must be one of the following strings "lrls", "krls", "krlsrf", "gpr"');
        else
            error('The option "algorithm" must be one of the following strings "lrls", "krls", "krlsrf", "gpr"');
        end
    end
    
    if ~isprop(opt, 'algorithm') || ~ok
        switch 1
            case regexp(opt.datatype, '\w*')
                algorithm = 'krls';
        end
        if ~ok
            fprintf('The selected algorithm is "%s".\n',algorithm);
        end
        opt.newprop('algorithm', algorithm);
    end
end

function filter = selectFilter(opt)
    if ~isprop(opt, 'filter')
        filter = 'tikh';
    elseif ischar(opt.filter)
        filter = unknelem('filter', opt.filter, '(tikh|land|nu|conjgrad|randtikh)', 'tikh');
    elseif isa(opt.filter, 'function_handle')
%         opt.newprop('filterfun', opt.filter);
%         opt.setting.filterfun = opt.filterfun;
%         filter = 'user-supplied';
        error('The option "filter" must be one of the following strings "tikh", "land", "nu", "conjgrad", "randtikh"');
    else
        error('The option "filter" must be one of the following strings "tikh", "land", "nu", "conjgrad", "randtikh"');
    end
    opt.newprop('filter', filter);
end

function chooseMatches(opt)
    if isequal(opt.datatype, 'vector')
        algorithm = selectAlgo(opt);
        if isequal(opt.algorithm, 'lrls')
            kernel = selectKernel(opt, true);
        else
            kernel = selectKernel(opt, false);
        end
    else
        kernel = selectKernel(opt, false);
        algorithm = selectAlgo(opt);
    end
    filter = selectFilter(opt);    
        
    type = 'primal';
    if opt.n <= opt.d
        type = 'dual';
    end
    algstr = [algorithm, '_', filter, '_', kernel, '_', type];
    switch 1
    case regexp(algstr, 'lrls_conjgrad_\w*_dual')
       algo = 'conjgraddual';
       pname = 'dual';
    case regexp(algstr, 'lrls_conjgrad_\w*_primal')
       algo = 'conjgradprimal';
       pname = 'primal';
    case regexp(algstr, 'lrls_pegasos_linear_\w*')
       algo = 'pegasos';
       pname = 'primal';
    case regexp(algstr, 'gpr_\w*_\w*_\w*')
       algo = 'gpregr';
       pname = 'gpregr';
       filter = allowelem('filter', 'the algorithm "gpr"', filter, 'tikh', 'tikh');
    case regexp(algstr, '(lrls_randtikh_linear_primal|krls_randtikh_linear_primal)')
       algo = 'primalr';
       pname = 'primal';
       algorithm = lrlsalg(algorithm);
    case regexp(algstr, '(krls_randtikh_\w*_\w*|lrls_randtikh_linear_dual)')
       algo = 'dualr';
       pname = 'dual';
       algorithm = krlsalg(algorithm);
    case regexp(algstr, '(lrls_land_linear_primal|krls_land_linear_primal)')
       algo = 'landweberprimal';
       pname = 'primal';
       algorithm = lrlsalg(algorithm);
    case regexp(algstr, '(krls_land_\w*_\w*|lrls_land_linear_dual)')
       algo = 'landweberdual';
       pname = 'dual';
       algorithm = krlsalg(algorithm);
    case regexp(algstr, '(lrls_nu_linear_primal|krls_nu_linear_primal)')
       algo = 'nuprimal';
       pname = 'primal';
       algorithm = lrlsalg(algorithm);
    case regexp(algstr, '(krls_nu_\w*_\w*|lrls_nu_linear_dual)')
       algo = 'nudual';
       pname = 'dual';
       algorithm = krlsalg(algorithm);
    case regexp(algstr, 'krlsrf_\w*_\w*_\w*')
       algo = 'randfeats';
       pname = 'randfeats';
       kernel = 'linear';
       filter = allowelem('filter', 'the algorithm "krlsrf"', filter, 'tikh', 'tikh');
    case regexp(algstr, '(lrls_\w*_\w*_primal|krls_\w*_linear_primal') % filter: tikhonov or not known 
       algo = 'primal';
       pname = 'primal';
       algorithm = lrlsalg(algorithm);
       sstr = sprintf('the algorithm "%s"', algorithm);
       filter = allowelem('filter', sstr, filter, '(tikh|randtikh|land|nu|user-supplied)', 'tikh');
    case regexp(algstr, 'user-supplied_\w*_\w*_\w*') % filter: tikhonov or not known
       algo = 'user'; 
       pname = 'dual';
       sstr = sprintf('the algorithm "%s"', algorithm);
       filter = allowelem('filter', sstr, filter, '(tikh|randtikh|land|nu|user-supplied)', 'tikh');
    case regexp(algstr, '(\w*_\w*_\w*_\w*|lrls_\w*_linear_dual)') % filter: tikhonov or not known
       algo = 'dual';
       pname = 'dual';
       if isequal(algorithm,'lrls')
           algorithm = krlsalg(algorithm);
       else
           algorithm = unknelem('algorithm', algorithm, 'krls', 'krls');
       end
       sstr = sprintf('the algorithm "%s"', algorithm);
       filter = allowelem('filter', sstr, filter, '(tikh|randtikh|land|nu|user-supplied)', 'tikh');
    end 
    
    opt.filter = filter;
    opt.setting.filter = filter;
    opt.algorithm = algorithm;
    opt.setting.algorithm = algorithm;
    opt.newprop('algname', algo);
    opt.setting.algname = algo;
    opt.newprop('predname', pname);
    opt.setting.predname = algo;
    opt.newprop('kernel', struct());
    if isequal(kernel, 'user-supplied')
        opt.kernel.func = opt.kernelfun;
    else
        opt.kernel.func = str2func(['kernel_', kernel]);
    end
    opt.kernelfun = kernel;
    opt.kernel.type = kernel;
    opt.setting.kernel = opt.kernel;
    
    fprintf('Final choice: "%s" algorithm, with filter "%s" and kernel "%s".\n', algorithm, filter, kernel);
end

function chooseCVAlg(opt)
    ok = false;
    
    if isprop(opt, 'partuning')
        if ischar(opt.partuning)
            ok = ~isequal('', unknelem('cross-validation method', opt.partuning, '(ho|loo)', ''));
        else
            error('Unknown value for variable "partuning"');
        end
    end
    
    if ~isprop(opt, 'partuning') || ~ok
        opt.newprop('partuning','ho');
    end
    
    fprintf('The selected cross-validation method is "%s".\n', opt.partuning);
    
    ok = false;
    if isprop(opt, 'perfm')
        if ischar(opt.perfm)
            if ~isequal('', unknelem('performance measure', opt.perfm, '(rmse|macroavg|precrec|gpregr)', ''))
                ok = true;
                sstr = [opt.problem, '_', opt.algorithm, '_', opt.perfm];
                switch 1
                    case regexp(sstr, '\w*_gpr_gpregr')
                    case regexp(sstr, '\w*_gpr_\w*')
                        warning('Have you considered to use the performance measure "gpregr" for the "gpr" algorithm?');
                    case regexp(sstr, '\w*_\w*_gpregr')
                        disp('The "gpregr" measure can be used only with the "gpr" algorithm.');
                        ok = false;
                    case regexp(sstr, 'regression_\w*_rmse')
                    case regexp(sstr, 'regression_\w*_\w*')
                        warning('Have you considered to use the performance measure "rmse" for this regression problem?');
                    case regexp(sstr, 'classification_\w*_(macroavg|precrec)')    
                    case regexp(sstr, 'classification_\w*_\w*')
                        warning('Have you considered to use the performance measure "macroavg" or "precrec" for this classification problem?');
                end
            end
        else
            error('Unknown value for variable "perfm"');
        end
    end
    
    if ~isprop(opt, 'perfm') || ~ok
        if isequal(opt.problem, 'regression')
            opt.newprop('perfm','rmse');
        else
            opt.newprop('perfm', 'macroavg');
        end
    end
    
    fprintf('The selected performance measure is "%s".\n', opt.perfm);
end

function pars = selectPars(opt)
    ok = false;
    pars = 'all';
    
    if isprop(opt, 'pars')
        if ischar(opt.pars)
            if ~isequal('', unknelem('"pars" option value', opt.pars, '(none|reg|ker|all)', ''))
                pars = opt.pars;
                switch opt.pars
                    case 'none'
                        ok = isprop(opt, 'kerpar') && isprop(opt, 'regpar');
                    case 'reg'
                        ok = isprop(opt, 'kerpar') || isequal(opt.datatype, 'kernel') || isequal(opt.kernelfun, 'linear') || isequal(opt.algorithm, 'lrls');
                        if isprop(opt, 'regpar')
                            warning('You are asking to find the regularization parameter by cross-validation, but you already specified it.');
                        end
                    case 'ker'
                        ok = isprop(opt, 'regpar');
                        if isprop(opt, 'kerpar')
                            warning('You are asking to find the kernel parameter by cross-validation, but you already specified it.');
                        end
                    case 'all'
                        if isprop(opt, 'kerpar') || isprop(opt, 'regpar')
                            warning('You are asking to find both kernel and regularization parameters by cross-validation, but you already specified one of them.');
                        end
                        ok = true;
                end
            end
        end
    end
   
    if ~isprop(opt, 'pars') || ~ok
        if isprop(opt, 'kerpar') && isprop(opt, 'regpar')
            pars = 'none';
        elseif ~isprop(opt, 'regpar') && (isprop(opt, 'kerpar') || isequal(opt.algorithm, 'lrls') || isequal(opt.kernelfun, 'linear') || isequal(opt.datatype, 'kernel'))
            pars = 'reg';
        elseif isprop(opt, 'regpar') && ~isprop(opt, 'kerpar')
            pars = 'ker';
        else
            pars = 'all';
        end
        opt.newprop('pars', pars);
    end
    switch pars
        case 'reg'
            disp('The regularization parameter will be chosen by cross-validation.');
        case 'ker'
            disp('The kernel parameter will be chosen by cross-validation.');
        case 'all'
            disp('Both the kernel and regularization parameters will be chosen by cross-validation.');
    end
end


function modelselection(opt)
    %algo = 'all', 'unkn', 'lrls', 'krls', 'gpr', 'krlsrf',  'krlsny', ...
    %filt = 'all', 'unkn', 'tikh', 'pcr', ...
    %kern = 'all', 'unkn', 'none', 'linear', 'polynomial', 'gaussian', 'user' ...
    %pars = 'all', 'unkn', 'reg', 'ker', ...
    %perfm = 'all', 'unkn', 'rmse', 'macroavg', 'logprob' ...
    %partuning = 'all', 'unkn', 'kcv', 'ho', 'loo'
    %split = 'all', 'unkn', 'none' ...
    %algo_filt_kern_pars_perfm_cvalg_split

    %%TODO ,fixsigmalambda,hodualr,hoprimalr,
    
    selectPars(opt);
    
    if isprop(opt, 'regrange')
        opt.newprop('paramsel.guesses', opt.regrange);
    end
    
    if isprop(opt, 'kerrange')
        opt.newprop('kernel.kerrange', opt.kerrange);
    end
    
   if isprop(opt, 'regpar')
        opt.newprop('paramsel.lambda', opt.regpar);
        opt.newprop('paramsel.lambdas', opt.regpar);
        opt.newprop('paramsel.guesses', opt.regpar);
    end
    
    if isprop(opt, 'kerpar')
        opt.newprop('kernel.kerrange', opt.kerpar);
        opt.newprop('paramsel.sigma', opt.kerpar);
    end
    
    if ischar(opt.perfm)
        opt.newprop('hoperf', str2func(['perf_', opt.perfm])); 
    else
        opt.newprop('hoperf', opt.perfm);
    end
    opt.newprop('paramsel.hoperf', opt.hoperf);
    opt.newprop('paramsel.nholdouts', opt.nholdouts);
    
    if ~isequal(opt.algname, 'user')
        opt.newprop('paramsel.optimizer', str2func(['rls_', opt.algname]));
    else
        opt.newprop('paramsel.optimizer', opt.algname);
    end
    
    %algo_filt_kern_pars_perfm_cvalg_split
    sstr = [opt.algorithm, '_', opt.filter, '_', opt.kernelfun, '_',...
        opt.pars, '_', opt.perfm, '_', opt.partuning];
    switch 1
        case regexp(sstr, '\w*_\w*_\w*_none_\w*_\w*')
            cvalgo = 'none';
        case regexp(sstr, 'lrls_pegasos_linear_reg_\w*_\w*')
            %subsize, calibfile, hoperf, singlelambda
            cvalgo = 'calibratesgd';
        case regexp(sstr, 'lrls_tikh_linear_reg_\w*_ho')
            % nlambda, smallnumber, hoperf, nholdouts
            cvalgo = 'hoprimal';
        case regexp(sstr, 'lrls_randtikh_linear_reg_\w*_ho')
            % nlambda, smallnumber, hoperf, nholdouts
            cvalgo = 'hoprimalr';
    	case regexp(sstr, 'lrls_tikh_linear_reg_\w*_loo')
            %nlambda, smallnumber
            cvalgo = 'loocvprimal';
    	case regexp(sstr, 'lrls_\w*_linear_reg_\w*_ho')
            if ~isprop(opt, 'paramsel') || ~isfield(opt.paramsel, 'guesses');
                error('A range for the regularization parameter is needed (option "regrange")');
            end
            cvalgo = 'bfprimal';
        case regexp(sstr, 'gpr_tikh_\w*_reg_\w*_loo')
            %nlambda, lambdamin, lambdamax
            cvalgo = 'loogpregr';
         case regexp(sstr, 'gpr_tikh_\w*_reg_gpregr_\w*')
             %lambdamin, lambdamax, nlambda
             cvalgo = 'gpregrLambdaGrid';
        case regexp(sstr, 'gpr_tikh_\w*_reg_\w*_ho')
            %nlambda, lambdamin, lambdamax
            cvalgo = 'hogpregr';
         case regexp(sstr, 'gpr_tikh_rbf_\w*_gpregr_\w*')
             %lambdamin, lambdamax, nlambda, sigmamin, sigmamax, nsigma, hoperf
             cvalgo = 'gpregrSigLamGrid';
        case regexp(sstr, 'gpr_tikh_rbf_\w*_\w*_ho')
            %nlambda, nsigma, lambdamin, lambdamax, sigmamin, sigmamax, hoperf
            cvalgo = 'siglamhogpregr';
        case regexp(sstr, 'gpr_tikh_rbf_\w*_\w*_loo')
            %nlambda, nsigma, lambdamin, lambdamax, sigmamin, sigmamax, hoperf
            cvalgo = 'siglamloogpregr';
        case regexp(sstr, 'krlsrf_tikh_linear_reg_\w*_ho')
            %nlambda, smallnumber, hoperf, nholdouts, randfeats.samplesize,
            %randfeats.D, 
            cvalgo = 'horandfeats';
        case regexp(sstr, 'krls_tikh_\w*_reg_\w*_loo')
            %kernel.type, nlambda, hoperf
            cvalgo = 'loocvdual';
        case regexp(sstr, 'krls_tikh_rbf_\w*_\w*_loo')
            %kernel.type, nlambda, hoperf, nsigma, sigmamin, sigmamax (+loocvdual)
            cvalgo = 'siglam';
         case regexp(sstr, 'krls_tikh_\w*_\w*_\w*_loo')
             %kernel.type, nlambda, hoperf, (+loocvdual)
             cvalgo = 'lootikhkrls';
        case regexp(sstr, 'krls_tikh_\w*_reg_\w*_ho')
            %nlambda, smallnumber, hoperf, nholdouts, kernel.type
            cvalgo = 'hodual';
        case regexp(sstr, 'krls_randtikh_\w*_reg_\w*_ho')
            %nlambda, smallnumber, hoperf, nholdouts, kernel.type
            cvalgo = 'hodualr';
        case regexp(sstr, 'krls_\w*_\w*_reg_\w*_ho')
            cvalgo = 'bfdual';
        case regexp(sstr, 'krls_tikh_rbf_\w*_\w*_ho')
            %kernel.type, nlambda, hoperf, nsigma, sigmamin, sigmamax (+hodual)
            cvalgo = 'siglamho';
        case regexp(sstr, 'krls_tikh_\w*_\w*_\w*_ho')
            cvalgo = 'hokrlstikh';
         case regexp(sstr, 'krls_\w*_\w*_\w*_\w*_\w*')
            if ~isprop(opt, 'paramsel') || ~isfield(opt.paramsel, 'guesses');
                error('A range for the regularization parameter is needed (option "regrange")');
            end
            cvalgo = 'hokrls';   
%        case regexp(sstr, '\w*_\w*_\w*_\w*_\w*_\w*')
%           %TODO
%           cvalgo = 'hoall';
%         case regexp(sstr, '\w*_\w*_\w*_\w*_\w*_loo')
%             %TODO
%             cvalgo = 'looall';
    end
    
    opt.newprop('cvalgo', cvalgo);  
end

function elem = unknelem(kind, elem, expr, newelem)
    if ~isequal(regexp(elem, expr), 1)
        if isequal(newelem, '')
            fprintf('The %s "%s" is not known.\n', kind, elem);
        else
            fprintf('The %s "%s" is not known, "%s" has been chosen.\n', kind, elem, newelem);
        end
        elem = newelem;
    end
end

function elem = allowelem(kind, where, elem, expr, newelem)
    if ~isequal(regexp(elem, expr), 1)
        if isequal(newelem, '')
            fprintf('The %s "%s" is not allowed for %s.\n', kind, elem, where);
        else
            fprintf('The %s "%s" is not allowed for %s, "%s" has been chosen.\n', kind, elem, where, newelem);
        end
        elem = newelem;
    end
end
    
function algorithm = lrlsalg(algorithm)
    if ~isequal(algorithm,'lrls')
           fprintf(['The primal formulation is more efficient for this problem.'...
            ' Thus the "lrls" algorithm has been chosen\n']);
           algorithm = 'lrls';
    end
end

function algorithm = krlsalg(algorithm)
    if ~isequal(algorithm,'krls')
           fprintf(['The dual formulation is more efficient for this problem.'...
            ' Thus the "krls" algorithm with linear kernel has been chosen\n']);
           algorithm = 'krls';
    end
end









