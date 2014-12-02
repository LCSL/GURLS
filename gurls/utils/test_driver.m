function [] = test_driver(task,dirname,opt)

    disp(dirname)

    if nargin<3; opt = struct(); end    

    opt.hoMOnline = 0;
    opt.nlambda = load('nlambda.txt');
    opt.nsigma = load('nsigma.txt');
    opt.nholdouts = load('nholdouts.txt');
    opt.hoproportion = load('hoproportion.txt');
    opt.smallnumber = load('smallnumber.txt');
    opt.eig_percentage = load('eig_percentage.txt');
    tmp = importdata('singlelambda.txt');
    opt.singlelambda = str2func(tmp{1});
    tmp = importdata('hoperf.txt');
    opt.hoperf = str2func(['perf_' tmp{1}]);

    allfiles = struct2cell(dir(dirname));
    allfiles = allfiles(1,:);

    taskroot = regexp(task,'_','split');
    task = str2func(task);
    
    filein = [dirname '/in.txt'];
    fid = fopen(filein,'wt');
    fclose(fid);
    fileout = [dirname '/out.txt'];
    fid = fopen(fileout,'wt');
    fclose(fid);
    
    switch taskroot{1}
        case 'split'
            X = load('Xtr.txt');
            y = load('ytr.txt');
            vout = task(X,y,opt);
            indices = [vout{1}.tr' vout{1}.va'] -1; 
            ntr = length(vout{1}.tr);
            vout = struct('indices',indices,'lasts',ntr);
            saveopt(vout,'split',dirname,fileout)
        case 'paramsel'
            X = load('Xtr.txt');
            y = load('ytr.txt');
            opt = addopt(opt,'split',allfiles,dirname,filein);
            if isfield(opt,'split')
                split.tr = opt.split.indices(1:opt.split.lasts) +1;
                split.va = opt.split.indices(opt.split.lasts+1:end) +1;
                opt.split = split;
            end
            opt = addopt(opt,'kernel',allfiles,dirname,filein);
            opt = checkopt(opt);
            vout = task(X,y,opt);
            saveopt(vout,'paramsel',dirname,fileout)
        case 'rls'
            X = load('Xtr.txt');
            y = load('ytr.txt');
            opt = addopt(opt,'paramsel',allfiles,dirname,filein);
            opt = addopt(opt,'kernel',allfiles,dirname,filein);
            vout = task(X,y,opt);
            saveopt(vout,'optimizer',dirname,fileout)
        case 'kernel'
            X = load('Xtr.txt');
            y = load('ytr.txt');
            opt = addopt(opt,'paramsel',allfiles,dirname,filein);  
            opt = checkopt(opt);
            vout = task(X,y,opt);
            saveopt(vout,'kernel',dirname,fileout)
        case 'predkernel'
            X = load('Xtr.txt');
            y = load('ytr.txt');
            opt = addopt(opt,'paramsel',allfiles,dirname,filein);        
            opt = addopt(opt,'rls',allfiles,dirname,filein);             
            opt = addopt(opt,'kernel',allfiles,dirname,filein); 
            opt = checkopt(opt);
            vout = task(X,y,opt);
            saveopt(vout,'predkernel',dirname,fileout)
        case 'pred'
            X = load('Xtr.txt');
            y = load('ytr.txt');
            opt = addopt(opt,'rls',allfiles,dirname,filein);        
            opt = addopt(opt,'predkernel',allfiles,dirname,filein);      
            opt = addopt(opt,'kernel',allfiles,dirname,filein);
            opt = checkopt(opt);
            vout = task(X,y,opt);
            saveopt(vout,'pred',dirname,fileout)
        case 'perf'
            X = load('Xtr.txt');
            y = load('ytr.txt');
            opt = addopt(opt,'pred',allfiles,dirname,filein);        
            vout = task(X,y,opt);
            saveopt(vout,'perf',dirname,fileout)
        case 'conf'
            X = load('Xtr.txt');
            y = load('ytr.txt');
            opt = addopt(opt,'pred',allfiles,dirname,filein);        
            vout = task(X,y,opt);
            saveopt(vout,'conf',dirname,fileout)
    end
       
end

function opt = addopt(opt,optfield,allfiles,dirname,typesfile)

    if strcmp(optfield,'rls'), 
        filesin = allfiles(strmatch('optimizer',allfiles));
    else
        filesin = allfiles(strmatch(optfield,allfiles));
    end
    if length(filesin)>0;
        if isfield(opt,optfield);
            optstruct = opt.(optfield);
        else
            optstruct = struct();
        end
        for i =1:length(filesin);
            varin = importdata([dirname '/' filesin{i}]);
            if iscell(varin); 
                varin = varin{1}; 
                vartype = 'string';
            elseif all(size(varin)==1)
                vartype = 'number';
            else
                vartype = 'matrix';
            end
            
            reg = regexp(filesin{i},'-','split');
            if length(reg)==1;
                optstruct = varin;
            else
                reg = regexp(reg{2},'\.','split');
                optstruct.(reg{1}) = varin;
            end
            
            fid = fopen(typesfile,'a');
            fprintf(fid,'%s %s\n',filesin{i},vartype);
            fclose(fid);

        end
        opt.(optfield) = optstruct;
    end
end

function [] = saveopt(vout,optfield,dirname,typesfile)

    if ~isstruct(vout);
        filename = [dirname '/' optfield '.txt'];
        savevar(vout,filename,typesfile)      
    else
        names = fieldnames(vout);
        for i = 1:length(names);
            var = vout.(names{i});
            if iscell(var)
%                 var
%                 for c = 1:length(var);
%                     filename = [dirname '/' optfield '-' names{i} '-' num2str(c) '.txt'];
                %    savevar(var{c},filename,typesfile);
%                 end

                filename = [dirname '/' optfield '-' names{i} '.txt'];                                
                
                for c = 1:length(var);
                    tmp(c,:) = reshape(var{c}, 1, numel(var{c}));
                end
                savevar(tmp,filename,typesfile);
                clear tmp;
            else
                if and(iscell(var),numel(var)==1);
                    var = var{1};
                end
                filename = [dirname '/' optfield '-' names{i} '.txt'];
                savevar(var,filename,typesfile);
            end
        end
    end
end

function []  = savevar(var,filename,typesfile)
    if isnumeric(var)
        if or(size(strfind(filename, 'split-indices'), 2) > 0, size(strfind(filename, 'split-lasts'), 2) > 0)
            var = cast(var, 'uint32');
            dlmwrite(filename,var,'delimiter',' ');
        else
            dlmwrite(filename,var,'delimiter',' ','precision','%1.11e');
        end
        if all(size(var)==1)
            vartype = 'number';
        else
            vartype = 'matrix';
        end
    else
        fid = fopen(filename,'wt');
        fprintf(fid,'%s',var);
        fclose(fid);
        vartype = 'string';

    end
    filename = regexp(filename,'/','split');  
    filename = filename{end};
    fid = fopen(typesfile,'a');
    fprintf(fid,'%s %s\n',filename,vartype);
    fclose(fid);
    
end

function opt = checkopt(opt)
% Checking compatibility, transforming to new GURLS opt
if ~isa(opt, 'GurlsOptions')
    warning('Compatibility mode with GURLS 1.0');
    opt = GurlsOptions(opt);
end
end
