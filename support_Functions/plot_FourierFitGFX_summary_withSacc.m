function plot_FourierFitGFX_summary_withSacc(cfg, dataIN,GFX_FourierNull)
%%
% This function plots the poster / slide summary, showing binned GFX
% with overlayed Fourier Series model (best fit).
% all 3 DVs are side by side

% called from j9__Plot_Data_FourierFits_GFX


GFX_headY = cfg.HeadData;

gaitfield = {'gc', 'doubgc'};
gaitprint = {'gait', 'stride'};
binfield = {'','_binned'};
nGaits_toPlot= cfg.nGaits_toPlot;
usebin=cfg.usebin;
ppantData=[];
plotHead=[];
% 3 colours 1 per DV
barCols= {'b', 'r', 'm', 'k'};
barCols= {'k','k'};
fntsize=12;

    iLR=3;
%% wrangle data to plot:

%which field of the datastructure to plot?
if strcmp(cfg.DV, 'RT')
    usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_rts'];
    ylabis = 'RT';
elseif strcmp(cfg.DV, 'Accuracy')
    usefield = [gaitfield{nGaits_toPlot} '_binned_Acc'];
    if ~cfg.norm
        ylabis=  cfg.DV;
    else
        ylabis = [cfg.DV 'norm: ' cfg.normtype];
    end
elseif strcmp(cfg.DV, 'Counts')
    usefield = [gaitfield{nGaits_toPlot} '_binned_counts'];
    ylabis = [cfg.type ' ' cfg.DV];

    if strcmp(cfg.type, 'Saccade')

        usefield = [gaitfield{nGaits_toPlot} '_sacc_all_binned_counts'];
        
    end
end

%collate data across subjs:
for isub= 1:size(dataIN,1)

    ppantData(isub,:)= dataIN(isub,iLR).(usefield);

%     if nGaits_toPlot==2
%         plotHead(isub,:) = GFX_headY(isub).([gaitfield{nGaits_toPlot} '_LRL']);
%     else
        plotHead(isub,:) = GFX_headY(isub).(gaitfield{nGaits_toPlot});
%     end
end
%% if normON , normalize as appropriate

if cfg.norm==1
    pM = nanmean(ppantData,2);
    meanVals= repmat(pM, 1, size(ppantData,2));


    if strcmp(cfg.normtype, 'absolute')
        data = ppantData - meanVals;
    elseif strcmp(cfg.normtype, 'relative')
        data = ppantData  ./ meanVals;
        data=data-1;
    elseif strcmp(cfg.normtype, 'relchange')
        data = (ppantData  - meanVals) ./ meanVals;
    elseif strcmp(cfg.normtype, 'normchange')
        data = (ppantData  - meanVals) ./ (ppantData + meanVals);
    elseif strcmp(cfg.normtype, 'db')
        data = 10*log10(ppantData  ./ meanVals);
    end

    ppantData= data;
end


%% other specs:
if nGaits_toPlot==1

    pidx= cfg.pidx1;
    ftnames= {'LR', 'RL', 'combined'};
else
    pidx= cfg.pidx2;
    ftnames= {'LRL', 'RLR', 'combined'};
end

%note that pidx is adjusted to all datapoints, if not using the  binned version.
if usebin==0
    pidx=1:size(ppantData,2);
end

% set xticks at approx centre points of the bins.
mdiff = round(mean(diff(pidx)./2));
xvec = pidx(1:end-1) + mdiff;

%% extracted fourier fits per shuffled series (calculated in j8_' loaded above)

%Hzspace% loaded above
fits_Rsquared_obsrvd=GFX_FourierNull.([cfg.type 'Ons_' usefield '_fitsRsq_Obs']);
fits_Rsquared_shuffCV=GFX_FourierNull.([cfg.type 'Ons_' usefield '_fitsRsq_ShuffCV']);

figure(1);
if cfg.iDV==1
    clf;
    set(gcf, 'units', 'normalized', 'position', [.05 .05 .9 .5], 'color', 'w');
end
% subplot(1,5, cfg.iDV)
subplot(1,3, cfg.iDV)

% gM = squeeze(mean(ppantData));
gM = squeeze(nanmean(ppantData));
stE = CousineauSEM(ppantData);
hold on;
% finely sampled bar, each gait "%" point.
bh=bar(xvec, gM);
hold on;
errorbar(xvec, gM, stE, ...
    'color', 'k',...
    'linestyle', 'none',...
    'linew', 2);
bh.FaceColor = [.9 .9 .9];
bh.EdgeColor = barCols{cfg.iDV};

sdrange = max(gM) - min(gM);
ylim([min(gM)-.5*sdrange max(gM)+1*sdrange])
ylsat = get(gca, 'ylim');


midp=xvec(ceil(length(xvec)/2));

%% apply best fit (overlay)
[f,gof]= fit(xvec',gM',  'fourier1');

% plot:
hold on;
%             yyaxis right
if cfg.iDV>1
h=plot(f, xvec, gM);%,
h(2).LineWidth = 2;
h(2).Color = barCols{cfg.iDV};
%%
%treat max xvec as our full 'period'
fitperiod = f.w;
%convert to period per samples.
% include period and Rsquared
%treat max xvec as our full 'period'
Hzapp = xvec(end)/ (2*pi/(f.w));
legdetails = [sprintf('%.2f', Hzapp) ' cp' gaitprint{nGaits_toPlot}(1) ', R^2 = ' sprintf('%.2f', gof.rsquare) ];
legend(h(2), legdetails, 'fontsize', fntsize-2, 'autoupdate', 'off', 'Location', 'South')
end
ylabel(ylabis)
set(gca,'fontsize', fntsize, 'xtick', [1, midp, xvec(end)], 'XTickLabels', {'0', '50', '100%'})
xlabel([ cfg.type ' onset (% ' gaitprint{nGaits_toPlot} '-cycle )']);%


%%
% overlay head position?
% hold on
% yyaxis right
% stEH= CousineauSEM(plotHead);
% pH= nanmean(plotHead,1);
% sh= shadedErrorBar(1:size(plotHead,2), pH, stEH,'k',1) ;
% sh.mainLine.LineStyle='-';
% sh.mainLine.LineWidth = 4;
% legp = sh.mainLine;
% set(gca,'ytick', [], 'fontsize',fntsize, 'YColor', 'k');
% 
% ylim([-.2 .05])
%>>>>> END DV with fits

%%  also show fit strength across frequencies: 
% subplot(1,5,5);
subplot(1,3,3);

hold on;

pO=plot(cfg.Hzspace(1:length(fits_Rsquared_obsrvd)), fits_Rsquared_obsrvd,'color', barCols{cfg.iDV}, 'linew', 3);

ylabel('R^2');
xlabel('Frequency (cycles-per-stride)')


legHZ(cfg.iDV) = pO;
hold on;
ph=plot(cfg.Hzspace(1:length(fits_Rsquared_obsrvd)), fits_Rsquared_shuffCV(3,:), ':', 'linew', 2, 'color', barCols{cfg.iDV});
set(gca,'fontsize', fntsize)
% title('(forced) Fits per frequency', 'fontsize', 18)
xlim([.5 6])
ylim([0 .75]);



%%
% legend([legHZ(2), legHZ(1), legHZ(3), ph], {'Accuracy', 'RT', 'Response', '95% CI'},'fontsize', fntsize-2)
ylim([0 .9])

end % function
%%%% >>>>>>>>>>>>>>>>>>>>>>>

