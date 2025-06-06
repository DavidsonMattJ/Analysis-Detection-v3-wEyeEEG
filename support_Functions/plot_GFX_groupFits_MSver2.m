function plot_GFX_groupFits_MSver2(cfg,testData, GFX_FourierNull)
% % % helper function to plot a MS ready summary of the group fits to
% detection performance. 

% will include an example trial, and target onset coutns over the gait.

% !! only the gait (linked steps),in this version:

% called from j9__Plot_Data_FourierFits_GFX

% % % % set some parameters:

GFX_headY = cfg.HeadData;

gaitfield = {'gc', 'doubgc'};
gaitprint = {'step', 'stride'};
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
for igait= 2%1:2 % both gait versions:

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
subplot(3,4,3)
hold on
% yyaxis right
stEH= CousineauSEM(plotHead);
pH= nanmean(plotHead,1);

%resize for plot?
pHs= imresize(pH, [1,100]);
stEHs= imresize(stEH, [1,100]);
sh= shadedErrorBar(1:100, pHs, stEHs,{'color','k','linew',1}) ;
sh.mainLine.LineStyle='-';
sh.mainLine.LineWidth = 4;
legp = sh.mainLine;
set(gca,'fontsize',cfg.fntsize, 'YColor', 'k');

% ylim([-.2 .05]) % if detrended
ylabel('detrended Head height')
xlabel( '% stride-cycle');

%% add swing/stance phases?
 markersAt = ceil([15, 75 115 175]/2); % note that 40-60 is the usual stance-swing split.
            yat = get(gca, 'ylim');
            %make room for data:
            ylim(yat)
            ylim([yat(1)-.005 yat(2)])
            Yverts = [yat(1) - .005, yat(1)];
            Xverts = [ 0, markersAt, 100]; 
            pos1 = [0 Yverts(1) markersAt(1) .005]; % x y width height
            pos2 = [markersAt(1) Yverts(1) markersAt(2)-markersAt(1) .005];
           pos3 = [markersAt(2) Yverts(1) markersAt(3)-markersAt(2), .005];
           pos4 = [markersAt(3) Yverts(1) markersAt(4)-markersAt(3), .005];
           pos5 = [markersAt(4) Yverts(1) 100-markersAt(4), .005];


            h=rectangle('Position', pos1, 'FaceColor',[.8 .8 .8]);
           h=rectangle('Position', pos2, 'FaceColor','w');
           h=rectangle('Position', pos3, 'FaceColor',[.8 .8 .8]);
           h=rectangle('Position', pos4, 'FaceColor','w');
           h=rectangle('Position', pos5, 'FaceColor',[.8 .8 .8]);
          
           text(mean(markersAt(1:2)), Yverts(1), 'Swing', 'VerticalAlignment','bottom', 'HorizontalAlignment','center', 'FontSize',8);          
           text(mean(markersAt(2:3)), Yverts(1), 'Stance', 'VerticalAlignment','bottom', 'HorizontalAlignment','center', 'FontSize',8);
          
           text(mean(markersAt(3:4)), Yverts(1), 'Swing', 'VerticalAlignment','bottom', 'HorizontalAlignment','center', 'FontSize',8);          
%            


%% 
plot([markersAt(2) markersAt(2)], [-.02 -.01], 'k--', 'linew',1);
text(markersAt(2)+3, -.005, 'approx.', 'HorizontalAlignment','right')
text(markersAt(2)+3, -.0075, 'heel-strike', 'HorizontalAlignment','right')

plot([markersAt(3) markersAt(3)], [-.02 -.01], 'k--', 'linew',1);
text(markersAt(3)+3, -.005, 'approx.', 'HorizontalAlignment','left')
text(markersAt(3)+3, -.0075, 'toe-off', 'HorizontalAlignment','left')
%%

% set xticks at approx centre points of the bins.
mdiff = round(mean(diff(pidx)./2));
xvec = pidx(1:end-1) + mdiff; % xaxis vector
midp=xvec(ceil(length(xvec)/2)); % mid point:

subplot(3,4,4);
params=[];
params.xvec= xvec;
params.col ='k';
params.ScatCol = 'k';
params.scatSize=1;
[bh, Fh]=bar_wFourier(ppantData,params);
% [bh, Fh]=bar_wFourier(xvec, ppantData,'k');

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
% ylim([min(gM)-.5*sdrange max(gM)+1*sdrange])
ylim([0 20]);

%% apply best fit (overlay)
% [f,gof]= fit(xvec',gM',  'fourier1');
% or apply linear best fit. 
%
[f,gof]= fit(xvec',gM',  'a*x+b');

% plot:
hold on;
%             yyaxis right
h=plot(f, xvec, gM);%,
h(2).LineWidth = 2;
h(2).Color = [.2 .2 .2];
%
%treat max xvec as our full 'period'
% fitperiod = f.w;
%convert to period per samples.
% include period and Rsquared
%treat max xvec as our full 'period'
% Hzapp = xvec(end)/ (2*pi/(f.w));
% legdetails = [sprintf('%.2f', Hzapp) '(cps), R^2 = ' sprintf('%.2f', gof.rsquare) ];

% legdetails = [' R^2 = ' sprintf('%.2f', gof.rsquare) ];

legend off;
% legend(h(2), legdetails, 'fontsize', cfg.fntsize-2, 'autoupdate', 'off', 'Location', 'NorthEast')
end

xlabel( '% stride-cycle');
ylabel('Target Counts');
shg









%% now cycle through data types and plot each time:

subspots = [5,9;6,10;7,11];

ic=1;
for igait=2%1:2

    pidx = xindexes{igait};
    iLR=3;
    nGaits_toPlot= igait;

for id= 1:length(testData)


dataIN= testData{id};


%% wrangle data to plot:

%which field of the datastructure to plot?
if strcmp(cfg.DV(id), 'RT')
    usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_rts'];
    ylabis = 'Reaction time [sec]';
elseif strcmp(cfg.DV(id), 'Accuracy')
    usefield = [gaitfield{nGaits_toPlot} '_binned_Acc'];
    if ~cfg.norm
        ylabis=  cfg.DV(id);
        ylabis = 'Hit Rate';
    else
        ylabis = [cfg.DV(id) 'norm: ' cfg.normtype];
    end
elseif strcmp(cfg.DV(id), 'Counts')
    usefield = [gaitfield{nGaits_toPlot} '_binned_counts'];
    ylabis = [cfg.type(id) ' ' cfg.DV(id)];
    ylabis = 'Response counts';
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



subplot(3,4, subspots(ic,:))
params=[];
params.xvec= xvec;
params.col = barCols{id};
params.ScatCol = barCols{id};
params.scatSize=1;
[bh, Fh]=bar_wFourier(ppantData,params);

gM = squeeze(mean(ppantData));
stE = CousineauSEM(ppantData);
hold on;
% finely sampled bar, each gait "%" point.
bh=bar(xvec, gM);
bh.FaceColor = [.9 .9 .9];
bh.EdgeColor = barCols{id};
bh.EdgeAlpha=.5;
hold on;
errorbar(xvec, gM, stE, ...
    'color', [.5 .5 .5],...
    'linestyle', 'none',...
    'linew', .5);

sdrange = max(gM) - min(gM);
ylim([min(gM)-1*sdrange max(gM)+1*sdrange])
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
legdetails = [sprintf('%.2f', Hzapp) ' (cps), R^2 = ' sprintf('%.2f', gof.rsquare) ];

% legdetails = [' R^2 = ' sprintf('%.2f', gof.rsquare) ];


legend(h(2), legdetails, 'fontsize', cfg.fntsize-2, 'autoupdate', 'off', 'Location', 'NorthEast')

ylabel(ylabis)
set(gca,'fontsize', cfg.fntsize, 'xtick', [1, midp, xvec(end)], 'XTickLabels', {'0', '50', '100%'})


% adjust y ticks for RT:
if id==2
set(gca,'fontsize', cfg.fntsize, 'ytick', [.34:.01:.4]);

end
xlabel([ cfg.type{id} ' onset (% ' gaitprint{nGaits_toPlot} '-cycle )']);%
ic=ic+1;




%Hzspace% loaded above
fits_Rsquared_obsrvd=GFX_FourierNull.([cfg.type{id} 'Ons_' usefield '_fitsRsq_Obs']);
% fits_Rsquared_shuffCV=GFX_FourierNull.([cfg.type{id} 'Ons_' usefield '_fitsRsq_ShuffCV']);
fits_Rsquared_shuffCV=GFX_FourierNull.([cfg.type{id} 'Ons_' usefield '_fitsRsq_ShuffCV_max']);

subplot(3,4,[8,12]);
hold on;

% if the first, plot a grey version and withold handle for the legend
if id==1
%     pG= plot(cfg.Hzspace(1:length(fits_Rsquared_obsrvd)), fits_Rsquared_shuffCV(3,:),':','linew',2,'color', [.3 .3 .3]);
    pG= plot(xlim, [fits_Rsquared_shuffCV(3),fits_Rsquared_shuffCV(3)],':','linew',2,'color', [.3 .3 .3]);
    
    % also add target onset data:
    fits_Rsquared_obsrvd_T=GFX_FourierNull.(['TargetOns_doubgc_binned_counts_fitsRsq_Obs']);
    fits_Rsquared_shuffCV_T=GFX_FourierNull.(['TargetOns_doubgc_binned_counts_fitsRsq_ShuffCV']);

    pT= plot(cfg.Hzspace(1:length(fits_Rsquared_obsrvd)), fits_Rsquared_obsrvd_T,'color', [.2 .2 .2], 'linew', 3);
    
    ph=plot(cfg.Hzspace(1:length(fits_Rsquared_obsrvd)), fits_Rsquared_shuffCV_T(3,:), ':', 'linew', 2, 'color', [.2 .2 .2]);

end
    pO=plot(cfg.Hzspace(1:length(fits_Rsquared_obsrvd)), fits_Rsquared_obsrvd,'color', barCols{id}, 'linew', 3);
ylabel('R^2');
xlabel(['Frequency (cycles-per-' gaitprint{igait} ')'])


legHZ(id) = pO;
hold on;
% ph=plot(cfg.Hzspace(1:length(fits_Rsquared_obsrvd)), fits_Rsquared_shuffCV(3,:), ':', 'linew', 2, 'color', barCols{id});
    ph= plot(xlim, [fits_Rsquared_shuffCV(3),fits_Rsquared_shuffCV(3)],':','linew',2,'color',  barCols{id});

set(gca,'fontsize', cfg.fntsize)
% title('(forced) Fits per frequency', 'fontsize', 18)
xlim([.5 10])
xlim([.1 5]);
ylim([0 .75]);
%% plot below patch?
% create the patch
% x=cfg.Hzspace(1:length(fits_Rsquared_obsrvd));
% y=fits_Rsquared_shuffCV(3,:);
% patch([x fliplr(x)], [zeros(size(x)) fliplr(y)], [.8 .8 .8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');

patch([0 0 cfg.Hzspace(end) cfg.Hzspace(end)], [0,fits_Rsquared_shuffCV(3),fits_Rsquared_shuffCV(3),0] ,[.8 .8 .8], 'FaceAlpha', 0.2, 'EdgeColor', 'none');



if id==3 || id==6
    legend(pG, '95% CI')
end


%% where was above the CI?
ytest= fits_Rsquared_obsrvd;
% tCV = fits_Rsquared_shuffCV(3,:);
tCV = fits_Rsquared_shuffCV(3);
sigAt = find(ytest>tCV);
disp(['GFX exceed significance at :']);
disp(['Hz:' string(cfg.Hzspace(sigAt))]);
disp(['Max R2:' num2str(max(fits_Rsquared_obsrvd(sigAt)))]);




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



end


% % 
function [bh, Fh]= bar_wFourier(ppantData,params)

% called to plot the nice bar charts with fourier overlayed. returns handle
% for bar (bh) and Fourier fit (Fh)

gM = squeeze(mean(ppantData));
stE = CousineauSEM(ppantData);
hold on;
% finely sampled bar, each gait "%" point.
bh=bar(params.xvec, gM,'BarWidth',1); % 1 and bars will touch.
hold on;
errorbar(params.xvec, gM, stE, ...
    'color', params.col,...
    'linestyle', 'none',...
    'linew', 1);
bh.FaceColor = [.9 .9 .9];
bh.EdgeColor = params.col;%barCols{id};

sdrange = max(gM) - min(gM);
ylim([min(gM)-.5*sdrange max(gM)+1*sdrange])
ylsat = get(gca, 'ylim');


%% apply best fit (overlay)
[f,gof]= fit(params.xvec',gM',  'fourier1');

% plot:
hold on;
%             yyaxis right
Fh=plot(f, params.xvec, gM);%,
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

% add scatter?
 hold on;
    % adjust scatter for within participant comparisons:
   [~, data_s] = CousineauSEM(ppantData);
    for id=1:size(data_s,2)
        useD = data_s(:,id);
        jitter = zscore(rand(length(useD), 1))/40;

        sc=scatter(params.xvec(id)+jitter, useD);
        sc.LineWidth = 2;
        sc.SizeData = params.scatSize;
        sc.MarkerFaceColor = params.ScatCol;
        sc.MarkerEdgeColor = params.ScatCol;
        sc.MarkerFaceAlpha = .4;
        sc.MarkerEdgeAlpha=.4;
        %
    end

end