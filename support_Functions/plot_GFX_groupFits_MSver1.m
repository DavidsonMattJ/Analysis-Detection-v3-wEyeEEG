function plot_GFX_groupFits_MSver1(cfg,testData, GFX_FourierNull)
% % % helper function to plot a MS ready summary of the group fits to
% detection performance. 

% will include an example trial, and target onset coutns over the gait.


% called from j9__Plot_Data_FourierFits_GFX

% % % % set some parameters:

GFX_headY = cfg.HeadData;

gaitfield = {'gc', 'doubgc'};
gaitprint = {'step', 'gait'};
binfield = {'','_binned'};
nGaits_toPlot= cfg.nGaits_toPlot;
usebin=cfg.usebin;
ppantData=[];
plotHead=[];
% 3 colours 1 per DV
barCols= {'b', 'r', 'm'};
% fntsize=12;

iLR=3;

xindexes = {cfg.pidx1, cfg.pidx2};

figure(1);
clf;
set(gcf, 'units', 'normalized', 'position', [.01 .01 .95 .95], 'color', 'w');

% first plot the example trial (function below).
subplot(3,4,1:2);

plotextrial(cfg);
% % now plot the target onsets for the next panels:
dataIN= testData{1}; % target onset info:
ylabis = ['Target counts'];
for igait= 1:2 % both gait versions:

    usefield = [gaitfield{igait} '_binned_counts'];
    pidx = xindexes{igait};
   
    %collate data across subjs:
    ppantData=[];
    plotHead=[]
    for isub= 1:size(dataIN,1)

        ppantData(isub,:)= dataIN(isub,iLR).(usefield);      
        plotHead(isub,:) = GFX_headY(isub).([gaitfield{igait}]);
       
    end

% set up plot:

% set xticks at approx centre points of the bins.
mdiff = round(mean(diff(pidx)./2));
xvec = pidx(1:end-1) + mdiff; % xaxis vector
midp=xvec(ceil(length(xvec)/2)); % mid point:

subplot(3,4,2+igait);
[bh, Fh]=bar_wFourier(xvec, ppantData);

% tidy with some specifics:
Fh(2).Color= 'k';
Fh(2).Visible= 'off'
ylabel(ylabis)
set(gca,'fontsize', cfg.fntsize, 'xtick', [1, midp, xvec(end)], 'XTickLabels', {'0', '50', '100%'})
xlabel(['Target onset (% ' gaitprint{igait} '-cycle )']);%
lg= get(gca, 'legend');
lg.Visible= 'off';
% ylim([0 10*igait])

gM= nanmean(ppantData,1);
sdrange = max(gM) - min(gM);
ylim([min(gM)-.5*sdrange max(gM)+1*sdrange])


% overlay head position?
hold on
yyaxis right
stEH= CousineauSEM(plotHead);
pH= nanmean(plotHead,1);

%resize for plot?
pHs= imresize(pH, [1,100]);
stEHs= imresize(stEH, [1,100]);
sh= shadedErrorBar(1:100, pHs, stEHs,{'color','k','linew',1}) ;
sh.mainLine.LineStyle='-';
sh.mainLine.LineWidth = 4;
legp = sh.mainLine;
set(gca,'ytick', [], 'fontsize',cfg.fntsize, 'YColor', 'k');

% ylim([-.2 .05]) % if detrended

end
legend(sh.mainLine, 'Head height')












%% now cycle through data types and plot each time:

subspots = [5,6,7,9,10,11];

ic=1;
for igait=1:2

    pidx = xindexes{igait};
    iLR=3;
    nGaits_toPlot= igait;

for id= 1:length(testData)


dataIN= testData{id};


%% wrangle data to plot:

%which field of the datastructure to plot?
if strcmp(cfg.DV(id), 'RT')
    usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_rts'];
    ylabis = 'RT';
elseif strcmp(cfg.DV(id), 'Accuracy')
    usefield = [gaitfield{nGaits_toPlot} '_binned_Acc'];
    if ~cfg.norm
        ylabis=  cfg.DV(id);
    else
        ylabis = [cfg.DV(id) 'norm: ' cfg.normtype];
    end
elseif strcmp(cfg.DV(id), 'Counts')
    usefield = [gaitfield{nGaits_toPlot} '_binned_counts'];
    ylabis = [cfg.type(id) ' ' cfg.DV(id)];
end

ppantData=[];
plotHead=[];
%collate data across subjs:
for isub= 1:size(dataIN,1)

    ppantData(isub,:)= dataIN(isub,iLR).(usefield);

    plotHead(isub,:) = GFX_headY(isub).([gaitfield{nGaits_toPlot}]);
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







subplot(3,4, subspots(ic))

[bh, Fh]=bar_wFourier(xvec, ppantData);

gM = squeeze(mean(ppantData));
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
bh.EdgeColor = barCols{id};

sdrange = max(gM) - min(gM);
ylim([min(gM)-.5*sdrange max(gM)+1*sdrange])
ylsat = get(gca, 'ylim');


midp=xvec(ceil(length(xvec)/2));

%% apply best fit (overlay)
[f,gof]= fit(xvec',gM',  'fourier1');

% plot:
hold on;
%             yyaxis right
h=plot(f, xvec, gM);%,
h(2).LineWidth = 2;
h(2).Color = barCols{id};
%%
%treat max xvec as our full 'period'
fitperiod = f.w;
%convert to period per samples.
% include period and Rsquared
%treat max xvec as our full 'period'
Hzapp = xvec(end)/ (2*pi/(f.w));
legdetails = [' R^2 = ' sprintf('%.2f', gof.rsquare) ];


legend(h(2), legdetails, 'fontsize', cfg.fntsize-2, 'autoupdate', 'off', 'Location', 'NorthEast')

ylabel(ylabis)
set(gca,'fontsize', cfg.fntsize, 'xtick', [1, midp, xvec(end)], 'XTickLabels', {'0', '50', '100%'})

xlabel([ cfg.type{id} ' onset (% ' gaitprint{nGaits_toPlot} '-cycle )']);%
ic=ic+1;




%Hzspace% loaded above
fits_Rsquared_obsrvd=GFX_FourierNull.([cfg.type{id} 'Ons_' usefield '_fitsRsq_Obs']);
fits_Rsquared_shuffCV=GFX_FourierNull.([cfg.type{id} 'Ons_' usefield '_fitsRsq_ShuffCV']);

subplot(3,4,8+(4*(nGaits_toPlot-1)));
hold on;

% if the first, plot a grey version and withold handle for the legend
if id==1
    pG= plot(cfg.Hzspace(1:length(fits_Rsquared_obsrvd)), fits_Rsquared_shuffCV(3,:),':','linew',2,'color', [.3 .3 .3]);
end
    pO=plot(cfg.Hzspace(1:length(fits_Rsquared_obsrvd)), fits_Rsquared_obsrvd,'color', barCols{id}, 'linew', 3);
ylabel('R^2');
xlabel(['Frequency (cycles-per-' gaitprint{igait} ')'])


legHZ(id) = pO;
hold on;
ph=plot(cfg.Hzspace(1:length(fits_Rsquared_obsrvd)), fits_Rsquared_shuffCV(3,:), ':', 'linew', 2, 'color', barCols{id});
set(gca,'fontsize', cfg.fntsize)
% title('(forced) Fits per frequency', 'fontsize', 18)
xlim([.5 10])
ylim([0 .9]);
%% plot below patch?
% create the patch
x=cfg.Hzspace(1:length(fits_Rsquared_obsrvd));
y=fits_Rsquared_shuffCV(3,:);
patch([x fliplr(x)], [zeros(size(x)) fliplr(y)], [.8 .8 .8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');





if id==3 || id==6
    legend(pG, '95% CI')
end







end % id
end % igait

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
%% extracted fourier fits per shuffled series (calculated in j8_' loaded above)

ylim([0 .9])











end
function plotextrial(cfg)

setmydirs_detectv3;
cd(procdatadir);
subjID= 'ABB';
itrial= 124;
cd(procdatadir)
lfile  = dir([pwd filesep subjID '*']);
load(lfile.name)
fntsize= cfg.fntsize;

%
plot(trialInfo(124).times, HeadPos(124).Y, 'k', 'LineWidth',2);
% axis tight;
xlabel('Trial time (sec)');
ylabel('Head height (m)');
set(gca,'fontsize', fntsize);
ylim([1.69 1.78])
%
% tOnsets = HeadPos(124).tr   
relT= find(trial_summaryTable.trial == itrial);
tOnsets = trial_summaryTable.targOnset(relT);
tResult = trial_summaryTable.targCor(relT);
hold on;
trialTimes = trialInfo(itrial).times;
 HeadData= HeadPos(itrial).Y;
txtHeight = 1.7;
 for itarg = 1:length(tOnsets)
   
     pAt = dsearchn(trialTimes, tOnsets(itarg)');

if tResult(itarg)==1% Hit.
%     ht=text(tOnsets(itarg), txtHeight, 'H', 'HorizontalAlignment', 'center')
    hl=plot(tOnsets(itarg), HeadData(pAt), 'ro', 'MarkerFaceColor', 'r');
    
    plot(tOnsets(itarg), 1.695, 'vk', 'MarkerFaceColor', 'k') 

else % Miss

%     mt=text(tOnsets(itarg), txtHeight, 'M', 'HorizontalAlignment','center')
    ml=plot(tOnsets(itarg), HeadData(pAt), 'ro', 'MarkerFaceColor', 'w');
    plot(tOnsets(itarg), 1.695, 'vk', 'MarkerFaceColor', 'w') 

end

% connect to base
plot([tOnsets(itarg) tOnsets(itarg)], [1.695 1.72], 'k-')
% plot(tOnsets(itarg), 1.7, 'dk', 'MarkerFaceColor', 'k') 
 end
 legend([hl,ml], {'Hit', 'Miss'})

 text(0.1, 1.71, {['Target'];['onset:']}, 'HorizontalAlignment', 'left', 'fontsize', fntsize-2);
 xlim([0 9])
 box off

end



% % 
function [bh, Fh]= bar_wFourier(xvec, ppantData)

% called to plot the nice bar charts with fourier overlayed. returns handle
% for bar (bh) and Fourier fit (Fh)

gM = squeeze(mean(ppantData));
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
% bh.EdgeColor = barCols{id};

sdrange = max(gM) - min(gM);
ylim([min(gM)-.5*sdrange max(gM)+1*sdrange])
ylsat = get(gca, 'ylim');


midp=xvec(ceil(length(xvec)/2));

%% apply best fit (overlay)
[f,gof]= fit(xvec',gM',  'fourier1');

% plot:
hold on;
%             yyaxis right
Fh=plot(f, xvec, gM);%,
Fh(2).LineWidth = 2;
%treat max xvec as our full 'period'
% fitperiod = f.w;
%convert to period per samples.
% include period and Rsquared
%treat max xvec as our full 'period'
% Hzapp = xvec(end)/ (2*pi/(f.w));
% legdetails = [sprintf('%.2f', Hzapp) ' cycles per, R^2 = ' sprintf('%.2f', gof.rsquare) ];
% legdetails = [' R^2 = ' sprintf('%.2f', gof.rsquare) ];

% legend(Fh(2), legdetails, 'fontsize', 8, 'autoupdate', 'off', 'Location', 'NorthEast')


end