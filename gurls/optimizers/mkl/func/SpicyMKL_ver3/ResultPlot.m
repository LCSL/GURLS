linstyles = {'-','-.','--'};
lnlen = length(linstyles);
markstyles = {'','o','>','^','<','d'};
mklen = length(markstyles);
colorstyles = {'b','m','r','c','k','g'};
collen = length(colorstyles);
font = 'Helvetica';
fontsize = 20;
linewidth = 3;
marksize = 13;

numsplit = 100;
numsplitx = 3;
end1 = floor(numsplit*0.6);
plotplace = [1:(end1*numsplitx)];
subplot(numsplit,numsplitx,plotplace);
hold off;        
for ii=1:size(DWEIGHT,2)
    ww = mod(ii,lnlen);
    ww2 = mod(floor(ii/lnlen),mklen);
    ww3 = mod(ii,collen);
    str = [linstyles{ww+1} markstyles{ww2+1} colorstyles{ww3+1}]; %'k'];
    plot(pam_array,DWEIGHT(:,ii),str,'LineWidth',linewidth,'MarkerSize',marksize);
    hold on;
end;
vv = kerneloptionvect{1};
clear LEG;
for ii=1:length(vv)
    LEG{ii} = sprintf('%g',vv(ii));
end;
legend(LEG,'FontName',font,'FontSize',fontsize);
if ismember(datatype,[1])
    xlabel('distance','FontName',font,'FontSize',fontsize)
elseif ismember(datatype,[2])
    xlabel('roughness','FontName',font,'FontSize',fontsize)
elseif ismember(datatype,[3])
    xlabel('frequency','FontName',font,'FontSize',fontsize)
end;
ylabel('kernel weight','FontName',font,'FontSize',fontsize)
set(gca,'FontName',font,'FontSize',fontsize*1.2);
xlim([min(pam_array) max(pam_array)]);


lenpam = length(pam_array);
if datatype == 1
    ylimvec = [-3 3];
elseif datatype == 2
    ylimvec = [-0.5 1];
elseif datatype == 3
    ylimvec = [0 1];
    xlimvec = [0 1.19];
end;
if ismember(datatype,[3])
    pamm = [min(pam_array) floor((min(pam_array)+max(pam_array))/2) max(pam_array)];
elseif ismember(datatype,[1 2]) 
    pamm = [min(pam_array) ((min(pam_array)+max(pam_array))/2) max(pam_array)];
end;
for ii=1:3
    subplot(numsplit,numsplitx,numsplitx*[(end1+15):numsplit-1]+ii);
    pam = pamm(ii); 

    param = setfield(param,pam_fld,pam);
    [xapp,yapp] = SyntheticDataGen(datatype,numsample,param);   

    hold off;
    plot(xapp((yapp>0),1),xapp((yapp>0),2),'bo','LineWidth',0.1,'MarkerSize',4); 
    hold on;
    plot(xapp((yapp<0),1),xapp((yapp<0),2),'r*','LineWidth',0.1,'MarkerSize',4); 
    ylim(ylimvec);
    if datatype == 3
        xlim(xlimvec);
    end;
    set(gca,'FontName',font,'FontSize',fontsize/1.5);
    grid on;
    if ismember(datatype,[1])
        str = sprintf('dist=%g',pam);
    elseif ismember(datatype,[2])
        str = sprintf('roughness=%0.3g',pam);
    elseif ismember(datatype,[3])
        str = sprintf('freq=%g',pam);
    end;
    xlabel(str,'FontName',font,'FontSize',fontsize/1.5);
end;
fname = sprintf('synthplot%d.eps',datatype);
print('-depsc',fname);