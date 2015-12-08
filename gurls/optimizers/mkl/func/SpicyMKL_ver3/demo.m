% Example of how to use the SpicyMKL for classification
%
%

clear all
close all
%addpath('./simpleMKL');

datacount = 0;

options.stopdualitygap = 1;
options.stopIneqViolation = 0;
options.tolOuter = 0.01;
options.tolInner = 0.0001;
options.loss = 'logit';


exps=[struct('datatype',1,...
            'numsample',1000,...
            'param', struct('noise',0.05,'freq',1.5,'dist',[]),...
            'kernelt',{{'gaussian'}},...
            'kerneloptionvect',{{[0.01 0.05 0.1 0.25 0.5 1]}},... 
            'pam_fld', 'dist',...
            'pam_array', linspace(10^(-2),1.2,20),...
            'C', 0.1),...
      struct('datatype',2,...
             'numsample',500,...
             'param',struct('noise',[]),...
             'kernelt',{{'gaussian'}},...
             'kerneloptionvect',{{[0.005 0.01 0.05 0.1 0.25 0.5 1]}},... 
             'pam_fld', 'noise',...
             'pam_array', logspace(-3,-0.7,30),...
             'C', 0.5),...
      struct('datatype',3,...
             'numsample',1000,...
             'param',struct('freq',[]),...
             'kernelt',{{'gaussian'}},...
             'kerneloptionvect' ,{{[0.01 0.025 0.05 0.1 0.25 0.5]}},... 
             'pam_fld', 'freq',...
             'pam_array', 2:2:22,...
             'C', 0.1)];

for jj=1:length(exps)
    datacount = datacount+1;
    pam_count = 0;
    variablevec={'all'};
    clear DWEIGHT;

    datatype  =exps(jj).datatype;
    numsample =exps(jj).numsample;
    pam_array =exps(jj).pam_array;
    pam_fld   =exps(jj).pam_fld;
    for pam = pam_array
        pam_count = pam_count + 1;

        param = setfield(exps(jj).param, pam_fld, pam);
        [xapp,yapp] = SyntheticDataGen(datatype,numsample,param);   

        [nbdata,dim]=size(xapp);
        
        kerneloptionvect=exps(jj).kerneloptionvect;
        [kernel,kerneloptionvec,variableveccell]=CreateKernelListWithVariable(variablevec,dim,exps(jj).kernelt,kerneloptionvect);
        [Weight,InfoKernel]=UnitTraceNormalization(xapp,kernel,kerneloptionvec,variableveccell);
        K=mklkernel(xapp,InfoKernel,Weight,options);
        [alpha,d,b,activeset,posind,params,story] = SpicyMKL(K,yapp,exps(jj).C,options);            
        disp(d)
        DWEIGHT(pam_count,:) = d;
    end;
%     figure(1);
%     hold off;        
%     plot(DWEIGHT,'LineWidth',3,'MarkerSize',13);
%     vv = kerneloptionvect{1};
%     clear LEG;
%     for ii=1:length(vv)
%         LEG{ii} = sprintf('%f',vv(ii));
%     end;
%     legend(LEG);

    figure,
    ResultPlot;

    fprintf('Press any key to continue.\n');
    pause;
end;



