function [y, performance] = gurls_test(model, Xtest, varargin)
    
    %% copy last process and modify to create a test process
    % default process to be used is always the last one in model.process 
    % if a single process exists, this falls back to original version 
    n_process = length(model.process);
    id_proc = n_process + 1;
    model.process{id_proc} = model.process{n_process};
    % ignore rls and paramsel jobs if these are part of seq
    model.process{id_proc}(strncmp(model.seq, 'paramsel:', 9)) = 3;
    model.process{id_proc}(strncmp(model.seq, 'rls:', 4)) = 3;
    model.process{id_proc}(strncmp(model.seq, 'kernel:', 7)) = 3;

    if numel(varargin) == 0
        model.process{end}(strncmp(model.seq, 'perf:', 5)) = 0;
        % model.process{end}(end) = 0;
        gurls(Xtest, [], model, id_proc);
        
        if nargout > 1
            error('In order to have the prediction error, the test labels must be provided');
        end
    else
        if numel(varargin) > 2
            error('Too many options');
        end
        
        ms3 = model.seq{end}; 
        if numel(varargin) == 2
            perfm = analyzePerfm(model, varargin{2});
            model.seq{end} = ['perf:' perfm];
        end
        if numel(varargin) >= 1
            ytrue = analyzeYtrue(model, varargin{1});
            model.process{end}(end) = 1;
            gurls(Xtest, ytrue, model, id_proc);
            
            switch model.seq{end}
                case 'perf:rmse'
                    performance = model.perf.rmse;
                case 'perf:macroavg'
                    performance = model.perf.acc;
                case 'perf:gpregr'
                    performance = model.perf.logprob;
                case 'perf:precrec'
                    performance = model.perf.ap;
                otherwise
                    performance = model.perf;
            end
            model.seq{end} = ms3;
            model.newprop('perf', struct());
        end
    end
    
    %% clean-up
    model.newprop('process', model.process(1:id_proc-1)); % remove dummy process entry?    
    y = convertFormat(model, model.pred);     
    % model.newprop('predkernel', struct()); % why do we need to delete the computed kernel?  
end

function ypred = convertFormat(model, y)

ypred = y;

if isstruct(model)
    
    if isfield(model.setting, 'yformat')
        switch model.setting.yformat
            case 'compact-twoclass'
                z = min(model.setting.labeldict)*ones(size(y,1),1);
                z(y>=0) = max(model.setting.labeldict);
                ypred = z;
            case 'compact-multiclass'
                [~, yp] = max(y,[], 2);
                ind = model.setting.labeldict;
                ypred = zeros(numel(yp),1);
                for i=1:size(ind)
                    ypred(yp==i) = ind(i);
                end
            case 'extended'
                % do nothing?
            case 'notinteger'
                % do nothing?
        end
    end
    
elseif isa(model, 'GurlsOptions')
    % do nothing
    % TO-DO: Id there any use for this in case model is a GurlsOpt object!
    % warning('GurlsOptions structure...nothing to convert');    
end
end

function ytrue = analyzeYtrue(model, y)

    T = model.setting.numoutputs;
    
    if isfield(model.setting, 'yformat')
        switch model.setting.yformat
            case 'compact-twoclass'
                T = 1;
                if size(y,2) ~= 1 || numel(unique(y)) > 2
                    error('The test labels are in a different format than the train labels.');
                end
                z = -ones(size(y,1),1);
                z(y == max(y)) = 1;
                ytrue = z;                
            case 'compact-multiclass'
                if size(y,2) ~= 1
                    error('The test labels are in a different format than the train labels.');
                elseif numel(setdiff(unique(y), model.setting.labeldict)) > 0
                    error('The test labels contain more classes than the train labels.');
                end
                ind = model.setting.labeldict;
                yp = zeros(numel(y),1);
                for i=1:size(ind)
                    yp(y==ind(i)) = i;
                end
                I = 2*eye(size(ind,1)) - 1;
                ytrue = I(yp,:);              
            case 'extended'
                ytrue = y;
            case 'notinteger'
                ytrue = y;
        end
    else
        ytrue = y;
    end
    
    if size(ytrue,2) ~= T
        error(['The test labels are %d-dimensional vectors, ', ...
            'while the training ones where %d-dimensional'], size(ytrue,2), T);
    end
end

function perfm = analyzePerfm(model, p)
    if ~isequal(regexp(p,'(rmse|macroavg|gpregr|precrec)'),1)
        error('Performance measure not recognized');
    else
        perfm = p;
    end
end
