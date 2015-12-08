function [alpha,d,b,activeset,supind,params,story] = SpicyMKL(K,yapp,C,options)
%
% a generalization of SpicyMKL, e.g. elastic net
%
% Usage:
%       [alpha,d,b,activeset,supind,params,story] = SpicyMKL(K,yapp,C,options)
%
% N: sample size,   M: number of kernels
%
% Input:
%    K        : N * N * M matrix. the (i,j,m)-element contains the
%               (i,j)-element of the m-th kernel gram matrix.
%    yapp     : N by 1 output signal vector. 
%    C        : regularization parameter. (large: strong, small: weak)
%               for l1 regularization C is a scalar: C|x|
%               for elasticnet, C is a two dimensional vector: C(1)|x| + C(2)x^2/2
%    options  : options which control the behavior of SpicyMKL
%       - loss: type of loss function: 'svm', 'logit', 'square'.
%                    (default: 'logit' for classification, 'square' for regression)
%                     svm: hinge loss for classification, max(0, 1-yf(x))
%                     logit: logistic regression, log(1+exp(- f(x)y))
%                     square: square loss, 0.5*(y - f(x))^2
%       - regname: type of regularization: 'l1', 'elasticnet'
%                    (default:'l1')                    
%       - outerMaxIter: maximum number of iteration of outer loop. (default 500)
%       - innerMaxIter: maximum number of iteration of inner loop. (default 500)
%       - stopdualitygap: 0 or 1. If 1, SpicyMKL employs duality gap for stopping criterion of outer loop.
%       - stopIneqViolation: 0 or 1. If 1, SpicyMKL employs violation of inequality for stopping criterion of outer loop.
%       - tolOuter: tollerance of stopping criteria of outer loop. (default 0.01)
%       - tolInner: tollerance of stopping criteria of inner loop. (default tolOuter/100)
%       - calpha: increment factor of gamma: gamma^(t+1)=calpha*gamma^(t). 
%                 (default 10)
%       - display: 1:no-display, 2:display outer loop, 3:display inner loop.
%                 (default 2)
%
% Output:
%     f(x) = \sum_{m \in activeset}(d(m)*km(x,:)*alpha) + b
%
%     alpha     :  N by 1 coefficient vector.
%     d         :  1 by M kernel weight vector.
%     b         :  bias.
%     activeset :  indices of kernels that are active ({m | d(m) is not zero}).
%     supind    :  indices of support vectors (when svm is used).
%     params    :  contains some variables.
%     story     :  contains history of primal objective, dual objective, 
%                     number of active kernels, and duality gap.
%                  - primalobj: primal objective
%                  - mod_dualobj: dual objective at modified dual variable
%                  - len_active: number of active kernels
%                  - dualitygap: duality gap
%
% (c) Taiji Suzuki and Ryota Tomioka: 
%     Department of Mathematical Informatics, The University of Tokyo, Japan. 
%     t_suzuki@mist.i.u-tokyo.ac.jp
%     tomioka@mist.i.u-tokyo.ac.jp


N = length(yapp);
M = size(K,3); 

if ~exist('options') 
    options = [];
end;
if ~isfield(options,'loss')
    if length(unique(yapp)) == 2
        options.loss = 'logit';
        yval = unique(yapp);
        yapp(yval==yval(1)) = -1;
        yapp(yval==yval(2)) = 1;
    else
        options.loss = 'square';
    end;
    fprintf('Please specify options.loss. We set options.loss=%s.\n',options.loss);
else
    bb = 0;
    loss_names = {'svm','logit','square'};
    for i = 1:length(loss_names)
        bb = bb | strcmp(options.loss,loss_names{i});
    end;
    if ~bb
        if length(unique(yapp)) == 2
            options.loss = 'logit';
            yval = unique(yapp);
            yapp(yval==yval(1)) = -1;
            yapp(yval==yval(2)) = 1;
        else
            options.loss = 'square';
        end;
        fprintf('Please specify options.loss correctly. We set options.loss=%s.\n',options.loss);
    end;    
end;
if ~isfield(options,'tolOuter')
    tol = 0.01;
else
    tol = options.tolOuter;
end;
if ~isfield(options,'tolInner')
    tolInner = tol/10000; 
else
    tolInner = options.tolInner;
end;
if ~isfield(options,'innerOptMeth')
    opt_method = 'Newton';
else
    opt_method = options.innerOptMeth;
end;
if ~isfield(options,'outerMaxIter')
    OuterMaxIter = 500;
else
    OuterMaxIter = options.outerMaxIter;    
end;
if ~isfield(options,'innerMaxIter')
    InnerMaxIter = 500;
else
    InnerMaxIter = options.InnerMaxIter;
end;
if ~isfield(options,'calpha')
    calpha = 10;
else
    calpha = options.calpha;
end;
if ~isfield(options,'stopIneqViolation')
    options.stopIneqViolation = 0;
end;
if ~isfield(options,'stopdualitygap')
    options.stopdualitygap = 1;
end;
if ~isfield(options,'display')
  options.display=2;
end
if ~isfield(options,'regname')
    options.regname = 'l1';
end;

lossname = options.loss;
regname = options.regname;


if strcmp(regname,'l1')
    regfunc = @(x,C)l1reg(x,C);
    regfuncdual = @(x,C)l1regdual(x,C);
    prox = @(x,C,eta)l1prox(x,C,eta);
    dprox = @(x,C,eta)l1dprox(x,C,eta);
elseif strcmp(regname,'elasticnet')
    if C(2) == 0
        regname = 'l1';
        C = C(1);
        regfunc = @(x,C)l1reg(x,C);
        regfuncdual = @(x,C)l1regdual(x,C);
        prox = @(x,C,eta)l1prox(x,C,eta);
        dprox = @(x,C,eta)l1dprox(x,C,eta);
    else
        regfunc = @(x,C)elasreg(x,C);
        regfuncdual = @(x,C)elasregdual(x,C);
        prox = @(x,C,eta)elasprox(x,C,eta);
        dprox = @(x,C,eta)elasdprox(x,C,eta);        
    end;
end;
    
story.primalobj = [];
story.mod_dualobj = [];
story.len_active = [];
story.dualitygap = [];
story.elapsedtime = [];
now_time = clock;
old_time = now_time;
elapsedtime = 0;
    
rho = -yapp/2;

mu = zeros(N,M);
cb = 0;

cgamma = ones(M,1)*10;
cgammab = 1;

if strcmp(regname,'elasticnet') && C(2)/C(1) > 5/95
    cgamma = ones(M,1)*1e+5;
    cgammab = 1e+5;
end;

if strcmp(options.loss,'svm')
    numt = 2;
    for i=1:numt
        ct{i} = ones(N,1)*10;  %\gamma_{\xi}
        clambda{i} = zeros(N,1);  %\xi
    end;
else
    clambda = [];
    ct = [];
end;

ck = Inf;

cbeta = 1/2;

u = rho*ones(1,M) + mu;
[v wj] = normKj(K,u);
qm = prox(wj.*cgamma,C,cgamma);

oneN = ones(N,1);
oneM = ones(1,M);

if strcmp(opt_method,'BFGS')
    Hk = eye(N);
end;

activeset = find(qm>0)'; %find(wj > C)';
supind = (1:N)';

     
%outer loop
for l=1:OuterMaxIter
    %inner loop
    for step = 1:InnerMaxIter    
        sumrho = sum(rho);
        if isempty(activeset)
            activeset = [];
        end;
        yrho = rho.*yapp;
        
        %compute gradient
        [fval,wvec] = funceval(lossname,wj,yapp,rho,yrho,sumrho,cgamma,cgammab,cb,qm,C,clambda,ct,regfunc);
        grad = gradient(lossname,yapp,yrho,rho,v,cgamma,activeset,wj,C,sumrho,qm,wvec,cgammab,cb);
            
        %compute Hessian
        if strcmp(opt_method,'Newton')
            dqm = dprox(wj.*cgamma,C,cgamma);
            alpha2 = (cgamma.*dqm.*wj - qm)./(wj.^3);
            alpha1 = qm./wj;
            if strcmp(lossname,'svm')
                Hessian = HessSVM(ct,yapp,wvec);    
                Hessian = HessAugMexbias(K,wj,activeset,v,C,cgamma,cgammab,Hessian,alpha1,alpha2); 
                if N < 2000
                    dk = - (Hessian + eye(N)*1e-8)\grad;                    
                else
                    dk = - pcg(Hessian + eye(N)*1e-8,grad,1e-4,200);
                end;
            elseif strcmp(lossname,'logit')
                Hessian = HessLogit(yrho); 
                Hessian = HessAugMexbias(K,wj,activeset,v,C,cgamma,cgammab,Hessian,alpha1,alpha2);    
                dk = - (Hessian)\grad;
                %dk = - pcg(Hessian,grad,1e-4,400); 
            elseif strcmp(lossname,'square')
                Hessian = HessSquare(length(rho));
                Hessian = HessAugMexbias(K,wj,activeset,v,C,cgamma,cgammab,Hessian,alpha1,alpha2);    
                dk = - (Hessian)\grad;
            end;
            graddotd = grad'*dk; 
        elseif strcmp(opt_method,'BFGS')
            dk = - Hk*grad;
            graddotd = grad'*dk; 
            if graddotd >=0 || step == 1
                dk = - grad;
                Hk = eye(N);
            else                    
                deltakBFGS = rho - old_rho_BFGS;
                gammakBFGS = grad - old_grad_BFGS;
                dgBFGS = (deltakBFGS'*gammakBFGS);
                Vk = eye(N) - deltakBFGS*(gammakBFGS/dgBFGS)';
                Hk = Vk*Hk*Vk' + deltakBFGS*(deltakBFGS/dgBFGS)';
            end;
            old_grad_BFGS = grad;    
            old_rho_BFGS = rho;
        end;
                    
        if graddotd > 0
            dk = - grad;
        end;
        
        step_size = 1; 
        old_fval = fval; 
        old_rho = rho;  
        old_wj = wj;        
        old_u = u; old_v = v;
        old_activeset = activeset; 
        old_yrho = yrho;
        
        rho = old_rho + step_size*dk;
        if  strcmp(lossname,'logit')
            yrho = rho.*yapp;
            if any(yrho <= -1 | yrho >= 0)
                yd = yapp.*dk;
                ss = min([(-1-old_yrho(yd<0))./yd(yd<0);-(old_yrho(yd>0)./(yd(yd>0)))])*0.99;
                step_size = min(step_size,ss);
                rho = old_rho + step_size*dk;
            end;
        end;
        u = rho*oneM + mu;
        [v wj] = normKj(K,u);
        qm = prox(wj.*cgamma,C,cgamma);
                        
        activeset = find(qm>0)'; 
        if isempty(activeset)
            activeset = [];
        end;
       
        %%%%%%Armijo%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sumrho = sum(rho);
        [fval] = funceval(lossname,wj,yapp,rho,yapp.*rho,sumrho,cgamma,cgammab,cb,qm,C,clambda,ct,regfunc);

        dir = rho - old_rho;

        tmp_activeset = union(activeset,old_activeset); 

        actdif = setdiff(tmp_activeset,activeset);
        old_u(:,actdif) = old_rho*ones(1,length(actdif))+mu(:,actdif);
        [old_v(:,actdif) old_wj(actdif)] = normKj(K,old_u,actdif);

        dirnorm = 0*wj;
        dirdotu = 0*wj;
        dirnorm(tmp_activeset) = dir'*(v(:,tmp_activeset)-old_v(:,tmp_activeset)); 
        dirdotu(tmp_activeset) = dir'*old_v(:,tmp_activeset);
        steplen = 1; 
        %search stepsize by Armijo's rule
        while fval > old_fval + steplen*0.1*(dir'*grad) && (steplen > 0)             
            steplen = steplen/2; 
            rho = old_rho + dir*steplen;
            wj = 0*wj;
            %max is for avoiding imagenary number
            wj(tmp_activeset) = sqrt(max(0,old_wj(tmp_activeset).^2 + 2*steplen*dirdotu(tmp_activeset) + steplen^2*(dirnorm(tmp_activeset)))); 
            qm = prox(wj.*cgamma,C,cgamma); %‚±‚±‚Í‘¬‚­‚·‚é—]’n‚ ‚è
            [fval] = funceval(lossname,wj,yapp,rho,yapp.*rho,sum(rho),cgamma,cgammab,cb,qm,C,clambda,ct,regfunc);                    
        end;
        %%%%%%end of Armijo%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if steplen ~= 1
            activeset = find(qm>0)'; 
            u = rho*oneM + mu;
            v(:,activeset) = normKj(K,u,activeset);
        end;

        if ~isreal(fval)
            norm(v)
        end;
        if options.display>=3
          fprintf('  [%d] fval:%g steplen:%g \n',step,fval,steplen);
        end
        if norm(old_rho - rho)/norm(old_rho) <= tolInner                        
            break;            
        end;
    end;    
    %end of inner loop (Newton or BFGS)
            
    sumrho = sum(rho);
    activeset = find(qm>0)'; %find(wj > C)';
    if isempty(activeset)
        activeset = [];
    end;
    yrho = yapp.*rho;
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%primal variable update
    mu = mu*0;
    if ~isempty(activeset)
        mu(:,activeset) = u(:,activeset).*(oneN*(qm(activeset)./wj(activeset)./cgamma(activeset))');
    end;    
    cb = cb + cgammab*sumrho;        

    %%computation of duality gap%%%%%%%%%%%%%%%%%%%%%% 
    alpha = -mu.*(oneN*cgamma');
    [vmu muwj] = normKj(K,alpha,activeset);
    
    modrho = rho;
    modrho = modrho - oneN*sum(modrho)/N;  
    [vv rhowj] = normKj(K,modrho*oneM,activeset);
    if ~isempty(activeset) 
        if strcmp(regname,'l1')
            modrho = modrho*min(1,C(1)/max(rhowj));  
            rhowj = rhowj*min(1,C(1)/max(rhowj));
        else 
            %modrho = rho;
        end;
    else
        %modrho = rho;
    end;
    aa = zeros(N,1);
    if strcmp(lossname,'svm')
        supind = find(clambda{2}==0);
    end;
    for j=activeset
        aa = aa - (K(:,supind,j)*alpha(supind,j));
    end;
    aa = aa + cb;
    if strcmp(lossname,'svm')
        aa = sum(max(aa.*yapp,-1) + 1);

        %%dual objective                  
        modyrho = yapp.*modrho;        
        mod_dualobj = -losssvm(modyrho) - sum(regfuncdual(rhowj,C));        
    elseif strcmp(lossname,'logit')
        aa = sum(log(1+exp(aa.*yapp))); 
        
        %%dual objective 
        modyrho = yapp.*modrho;
        modyrho = min(max(modyrho,-0.9999999),0.00000001);                
        mod_dualobj = - losslogit(modyrho)  - sum(regfuncdual(rhowj,C));        
    elseif strcmp(lossname,'square')
        aa = sum((aa + yapp).^2)*0.5;                 
        
        %%dual objective 
        mod_dualobj = -losssquare(yapp,modrho) - sum(regfuncdual(rhowj,C));        
    end;    
    %%primal objective 
    if ~isempty(activeset) 
        primalobj = aa + sum(regfunc(muwj,C));%sum(regfunc(muwj.*cgamma(activeset),C));        
    else
        primalobj = aa;
    end;
    
    story.primalobj = [story.primalobj primalobj];
    story.mod_dualobj = [story.mod_dualobj mod_dualobj];
    story.len_active = [story.len_active length(activeset)];
    dualitygap = (abs(primalobj - mod_dualobj)/abs(primalobj));
    story.dualitygap = [story.dualitygap dualitygap];
    now_time = clock;
    elapsedtime = elapsedtime+etime(now_time,old_time);
    story.elapsedtime = [story.elapsedtime elapsedtime];
    old_time = now_time;
    %%end of computation of duality gap%%%%%%%%%%%%%%%%%%%%%% 
    
    uu = u;
    if ~isempty(activeset) 
        uu(:,activeset) = u(:,activeset).*(oneN*min(1,C(1)./wj(activeset))');
    end;
    %hresid = (wj.*cgamma - qm).^2; 
    hresid = sqrt(sum((uu-rho*oneM).^2,1)); 
    maxgap = 0;
    if strcmp(lossname,'svm')
        work = abs(max( - 1 - yrho, - clambda{1}./ct{1}));
        ctI{1} = find( work > cbeta*ck );
        maxgap = max(max(work),maxgap); 
        work = abs(max( yrho, - clambda{2}./ct{2}));
        ctI{2} = find( work > cbeta*ck );
        maxgap = max(max(work),maxgap);        
    end;
    work = abs(sumrho);
    b_cb = ( work > cbeta*ck);
    maxgap = max(work,maxgap);

    I3 = find(hresid > cbeta*ck);
    maxgap = max(max(abs(hresid)),maxgap);

    if options.display>=2 && options.stopIneqViolation
      fprintf('[[%d]] maxgap:%g ck:%g',l,maxgap,ck);
    end
    if maxgap <= ck
        ck = maxgap;
    end;
    %%KKT stopping criterion%%%%%%%%%%%%%%%%%%%
    if ck < tol && options.stopIneqViolation
        break;
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    if strcmp(lossname,'svm')
        clambda{1} = max(clambda{1} + ct{1}.*(-yrho-1),0);
        clambda{2} = max(clambda{2} + ct{2}.*(yrho),0);
    end; 

    allind = 1:N;
    if strcmp(lossname,'svm')
        for i=1:2
            ct{i}(ctI{i}) = calpha * ct{i}(ctI{i});
            ctIcomp = setdiff(allind,ctI{i});
            ct{i}(ctIcomp) = calpha + ct{i}(ctIcomp);
        end;
    end;
    if b_cb
        cgammab = calpha * cgammab;
    else
        cgammab = calpha + cgammab;
    end;    
    cgamma(I3) = calpha * cgamma(I3);
    mu(:,I3) = mu(:,I3)/calpha;    
    cI3 = setdiff(1:M,I3);
    if ~isempty(cI3)
        cgamma(cI3) = calpha + cgamma(cI3);
        mu(:,cI3) = mu(:,cI3).*(oneN*((cgamma(cI3)-calpha)./cgamma(cI3))');
    end;
            
    u = rho*oneM + mu;    
    [v wj] = normKj(K,u);
    qm = prox(wj.*cgamma,C,cgamma);
    activeset = find(qm>0)';

    
    if options.display>=2
      if ~options.stopIneqViolation
          fprintf('[[%d]] ',l);
      end;
      fprintf('primal:%g dual:%g duality_gap:%g\n',primalobj,mod_dualobj,dualitygap);
    end
      
    if options.stopdualitygap && dualitygap < tol
        break;        
    end;       
end;
%end of outer loop

if options.display>=2
    fprintf('\n');
end;

params.cgamma = cgamma;
params.cgammab = cgammab;
params.wj = wj;
params.C = C;
params.loss = lossname;
if strcmp(lossname,'svm')
    params.ct = ct;
    params.clambda = clambda;
end;

work = zeros(size(qm));
work(activeset) = qm(activeset)./wj(activeset); 
d = work./(1-work./cgamma);
sumd = sum(d);
if sumd ~= 0
    d = d/sumd;
end;
alpha = -rho*sumd;

if strcmp(lossname,'svm')
    supind = find(clambda{2}==0);
else
    supind = (1:N)';
end;
b = - cb;


%common%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function grad = gradAug(v,cgamma,activeset,wj,C,grad,sumrho,cgammab,cb)       
    if ~isempty(activeset)
        for j = activeset   
            grad = grad + v(:,j)*(cgamma(j)*(wj(j)-C)/wj(j));
        end;
    end;
    if nargin >=7   
        grad = grad + (cgammab*sumrho+cb);
    end;

function [val,wvec] = funceval(lossname,wj,yapp,rho,yrho,sumrho,cgamma,cgammab,cb,qm,C,clambda,ct,regfunc)
    if strcmp(lossname,'svm')
        [val] = losssvm(yrho); 
        wvec{1} = max(0,clambda{1} + ct{1}.*(-1 - yrho));
        wvec{2} = max(0,clambda{2} + ct{2}.*(yrho));
    elseif strcmp(lossname,'logit')
        val = losslogit(yrho);
        wvec = [];
    elseif strcmp(lossname,'square')
        val = losssquare(yapp,rho);
        wvec = [];
    end;
    val = val - (sum(regfunc(qm,C)) + sum((qm.^2)./(2*cgamma)) - wj'*qm);
    val=  val + cgammab*sumrho^2/2 + cb*sumrho;
    if strcmp(lossname,'svm')
        val = val + sum((wvec{1}.^2 - clambda{1}.^2)./ct{1})/2 + sum((wvec{2}.^2 - clambda{2}.^2)./ct{2})/2;
    end;
    
function grad = gradient(lossname,yapp,yrho,rho,v,cgamma,activeset,wj,C,sumrho,qm,wvec,cgammab,cb) 
    if strcmp(lossname,'svm')
        grad = gradSVM(wvec,yapp);   
    elseif strcmp(lossname,'logit')
        grad = gradLogit(yrho,yapp);
    elseif strcmp(lossname,'square')
        grad = gradSquare(rho,yapp);      
    end;  
    if ~isempty(activeset)
        %grad = grad + sum(v(:,activeset).*(ones(size(v,1),length(activeset))*(qm(activeset)./wj(activeset))));
        for j = activeset   
            grad = grad + v(:,j)*(qm(j)/wj(j));
        end;
    end;
    %if nargin >=7   
        grad = grad + (cgammab*sumrho+cb);
    %end;
    
            
%svm%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function [HessBarrier] = HessSVM(ct,yapp,wvec)       
    HessBarrier = zeros(length(yapp)); 
    ind1 = (wvec{1}>0); ind2 = (wvec{2}>0);
    HessBarrier(ind1,ind1) = diag(ct{1}(ind1));
    HessBarrier(ind2,ind2) = HessBarrier(ind2,ind2) + diag(ct{2}(ind2));
    
function grad = gradSVM(wvec,yapp)       
    grad = ((-wvec{1} + wvec{2}).*yapp) + yapp;

function [val,wvec1,wvec2] = funcevalsvm(wj,yrho,sumrho,cgamma,cgammab,cb,clambda,ct,C)

    wvec1 = max(0,clambda{1} + ct{1}.*(-1 - yrho));
    wvec2 = max(0,clambda{2} + ct{2}.*(yrho));
    val = losssvm(yrho) + (cgamma'*(max(wj-C,0).^2))/2;
    val = val + sum((wvec1.^2 - clambda{1}.^2)./ct{1})/2 + sum((wvec2.^2 - clambda{2}.^2)./ct{2})/2;
    val = val + cgammab*sumrho^2/2 + cb*sumrho;

function val = losssvm(yrho)
    val = sum(yrho);
    
    
%logit%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function Hessian = HessLogit(yrho) 
    Hessian = diag(1./(-yrho.*(1+yrho)));    
    
function grad = gradLogit(yrho,yapp)       
    grad = yapp.*log((1+yrho)./(-yrho));
    
function [val] = funcevallogit(wj,yrho,cgamma,C, sumrho,cgammab,cb)
    if any(yrho>0 | yrho<-1)
      val = inf;
      return;
    end

    val = losslogit(yrho) + (cgamma'*(max(wj-C,0).^2))/2;
    
    if nargin >=5   
        val = val + cgammab*sumrho^2/2 + cb*sumrho;
    end;

function val = losslogit(yrho)
    val = sum((1+yrho).*log(1+yrho)-yrho.*log(-yrho));
        
%square%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function [val] = funcevalsquare(wj,yapp,rho,cgamma,C,sumrho,cgammab,cb)
    val = losssquare(yapp,rho) + (cgamma'*(max(wj-C,0).^2))/2;
    
    if nargin >=6   
        val = val + cgammab*sumrho^2/2 + cb*sumrho;
    end;

function val = losssquare(yapp,rho)
    val = 0.5*(rho'*rho + 2*rho'*yapp); 

function Hessian = HessSquare(N) 
    Hessian = eye(N);    
    
function grad = gradSquare(rho,yapp)       
    grad = rho + yapp;
    
    
%regularization term%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function reg = l1reg(x,C)
    reg = C*abs(x);    
function regd = l1regdual(x,C)
    regd = zeros(length(x),1);
    regd(x > C) = Inf;
function prx = l1prox(x,C,eta)
    prx = max(x - C*eta,0); 
function dprx = l1dprox(x,C,eta)    
    dprx = 1*(x > C*eta);

function reg = elasreg(x,C)
    reg = C(1)*abs(x) + (C(2)/2)*(x.^2);    
function regd = elasregdual(x,C)
    regd = (0.5/C(2))*(x.^2 - 2*C(1)*abs(x) + C(1)^2).*(x > C(1));
function prx = elasprox(x,C,eta)
    prx = max(x - C(1)*eta,0)./(1+C(2)*eta); 
function dprx = elasdprox(x,C,eta)    
    dprx = (1./(1+C(2)*eta)).*(x > C(1)*eta);
