


function  plot_PFX_FourierFits_overlayed(cfg, dataIN, PFX_FourierNull)
%% plot an overalyed summary of individual fits to a certain DV

 % called from j11__Plot_Data_FourierFits_PFX



 
figure(1); clf; set(gcf, 'color', 'w', 'units', 'normalized', 'position', [0.1 0.1 .75 .75]);
nsubs = length(cfg.subjIDs);


% plot_FourierFit(cfg,dataIN);
%%
GFX_headY = cfg.HeadData;

clf
plotn=1;
for nGaits_toPlot=2
pc=1; % plot counter
pspots = [1,2,3,4,5,6]; %suplot order
psubj= 'GFX'; % print ppid.
% both this and the next use the same figure function:
iLR=3;
gaitfield = {'gc', 'doubgc'};
gaitprint = {'gait', 'stride'};
binfield = {'','_binned'};
    
    legp=[]; % for legend
    ppantData=[];
    plotHead=[];
    
    %which field of datastructure to plot?
    if strcmp(cfg.DV, 'RT')
        usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_rts'];        
        ylabis = 'z(RT)';
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
    end
    
    %collate data:
    for isub= 1:size(dataIN,1)

        ppantData(isub,:)= dataIN(isub,iLR).(usefield);

        plotHead(isub,:) = GFX_headY(isub).(gaitfield{nGaits_toPlot} );
    end

    %% extract mean per subj, incase we want that ordering:
    meanDV=[];
    for isub=1:size(dataIN)
        meanDV(isub) = round(mean(dataIN(isub, iLR).([usefield])),2);
    end
     
    % if normON , normalize as appropriate
    
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
    
     

if cfg.reOrderbyMean==1 % reorder the PpantData in acending accuracy:
    
   [m,n] = sort(meanDV);
   ppantDatas= ppantData(n,:);
   ppantData= ppantDatas;
   meanAcc = m;
   subjIDs = subjIDs(n);
    PFX_FourierNull = PFX_FourierNull(1,n);
end
    
%% set up the Fourier model (restrict fit range to 0- 10 Hz).   
    % Declaring the type of fit.
    FitType = 'fourier1';
    % Creating and showing a table array to specify bounds.
    CoeffNames = coeffnames(fittype(FitType));    
    %set bounds for w
    CoeffBounds = array2table([-Inf(1,length(CoeffNames));...
        Inf(1,length(CoeffNames))],'RowNames',...
        ["lower bound", "upper bound"],'VariableNames',CoeffNames);

%specify range for fit, from 0 to 10 Hz:
   
        testw = 2*pi*10/100;
        
        CoeffBounds.w(1) = 0;
        CoeffBounds.w(2) = testw;
        
        %update fit opts settings
        
        %Update Fit Options setting.
        FitOpts = fitoptions('Method','NonlinearLeastSquares','Lower',table2array(CoeffBounds(1,:)),...
            'Upper',table2array(CoeffBounds(2,:)));



    %% other specs:
    if nGaits_toPlot==1
        
        pidx= cfg.pidx1;
        ftnames= {'LR', 'RL', 'combined'};
    else
        pidx= cfg.pidx2;
        ftnames= {'LRL', 'RLR', 'combined'};
    end
    
    %note that pidx is adjusted to all datapoints, if not using the  bin.
    if cfg.usebin==0
        pidx=1:size(ppantData,2);
    end
    
    %x axis:          %approx centre point of the binns.
            mdiff = round(mean(diff(pidx)./2));
            xvec = pidx(1:end-1) + mdiff;
    
    %% extracted fourier fits per shuffled series (calculated in j8_' loaded above)
    
    %Hzspace% loaded above
%      fits_Rsquared_obsrvd=GFX_FourierNull.([cfg.type 'Ons_' usefield '_' speedsAre{ispeed} '_fitsRsq_Obs']);
%      fits_Rsquared_shuffCV=GFX_FourierNull.([cfg.type 'Ons_' usefield '_' speedsAre{ispeed} '_fitsRsq_ShuffCV']);
    
    %% %%%%%%%%%%%%%%%%%%% first row of plot. 
    % Mean data, with errorbars head pos overlayed.
%     subplot(3,2,pspots(plotn))
%   
%     stEH= CousineauSEM(plotHead);
%     pH= nanmean(plotHead,1);
%    sh= shadedErrorBar(1:size(plotHead,2), pH, stEH,'k',1) ; 
%    sh.mainLine.LineStyle='-';
%    sh.mainLine.LineWidth = 4;
%     legp = sh.mainLine;
%     set(gca,'ytick', [], 'fontsize',cfg.fntsize-5);
%     lg=legend(legp, 'Head height [a.u.]', 'location', 'NorthEast', 'fontsize', cfg.fntsize-5);
%    %% lengthen yyaxis right so legend doesnt obscure
%     ytop = max(pH);
%     ylt= get(gca,'ylim');
%     sdrange = max(pH) - min(pH);
%    ylim([min(pH) 1.2*max(pH)])
%    
%     %%
%      midp=xvec(ceil(length(xvec)/2));
%     set(gca,'fontsize', cfg.fntsize, 'xtick', [1, midp, xvec(end)], 'XTickLabels', {'0', '50', '100%'})   
% %     xlabel([ cfg.type ' onset relative to ' gaitprint{nGaits_toPlot} '-cycle (%)']);%
%     
%     plotn= plotn+1;
%     
%     if nGaits_toPlot==1
%         continue % move on with plots.
%     end
    %% %%%%%%%%%%%%%%%%%%% second row of plot. 
    % Individual fourier fits (labelled per ppant).
    %% also prepare next overlay:
    
    % finely sampled bar, each gait "%" point.
   leg=[];
   allFits=[];
   allFreqs=[];
   for ippant = 1:size(dataIN,1)
       
       %restrict fit to shorter range(?)
       
    [f,gof]= fit(xvec',ppantData(ippant,:)',  'fourier1',FitOpts);
    
    % plot:
    hold on;
    subplot(2,2,1);   hold on;
    %             yyaxis right
    h=plot(f, xvec, ppantData(ippant,:)');%,
    h(1).Visible= 'off';
    h(2).Visible= 'off';
    % rescale.
    pfit = h(2).YData;
    pfitO = pfit-mean(pfit);%./max(abs(pfit)); % between 1 - 1.
%     pfit = (pfit./max(abs(pfit))) + ippant;
    pfitP = pfitO+ ippant*.2;
    
    subplot(2,2,[2,4]);  hold on; 
    pl=plot(h(2).XData, pfitP);
    leg(ippant) = pl;
    legend 'off';
    
    Hzapp = xvec(end)/ (2*pi/(f.w));
    
    % add initial, and mean Acc as label
    text(pidx(end), pfitP(end),  [cfg.subjIDs{ippant}(1:2) ' ' sprintf('%.2f', Hzapp)], 'fontsize', 15);
    
    subplot(2,2,1);  hold on; 
    plot(h(2).XData, pfitO); 
    legend off
    allFits(ippant,:) = pfitO;
    
    %%
    subplot(2,2,3); % plot correct order:
    tmpD= PFX_FourierNull(ippant).([cfg.type 'Ons_' usefield '_fitsRsq_Obs']);
    plot(cfg.Hzspace, tmpD, 'color',[.8 .8 .8])
    %%
        allFreqs(ippant,:) =tmpD;

   end
   %%
   subplot(2,2,1); %plot average of fits
   hold on; plot(h(2).XData, mean(allFits,1), 'linew', 8, 'color', 'k')
    %%
    [f,gof]= fit(h(2).XData',mean(allFits,1)',  'fourier1');
    pF=plot(f, h(2).XData, mean(allFits,1));
%     pF(2).LineWidth=8;
    %treat max xvec as our full 'period'
    fitperiod = f.w;
    %convert to period per samples.
    % include period and Rsquared
    %treat max xvec as our full 'period'
%     Hzapp = xvec(end)/ (2*pi/(f.w));
    Hzapp = h(2).XData(end)/ (2*pi/(f.w));
    
    title(['Best fit (to avg): ' sprintf('%.2f', Hzapp) ' cycles-per-' gaitprint{nGaits_toPlot} ' (cp' gaitprint{nGaits_toPlot}(1) ')']),
    legend off;
%     legdetails = [sprintf('%.2f', Hzapp) ' cp' gaitprint{nGaits_toPlot}(1) ', R^2 = ' sprintf('%.2f', gof.rsquare) ];
%     legend(h(2), legdetails, 'fontsize', fntsize, 'autoupdate', 'off')
    %%
    ylabel('Fit modulation');
    set(gca,'fontsize', cfg.fntsize, 'xtick', [1, 50, 100], 'XTickLabels', {'0', '50', '100%'})
    xlabel([ cfg.type ' onset as % of ' gaitprint{nGaits_toPlot} '-cycle ']);%
    %%
   
    
    subplot(2,2,3);
    plot(cfg.Hzspace, mean(allFreqs,1), 'linew', 2);
    ylabel('R^2');
    xlabel('cps')
    set(gca, 'fontsize', 15);
    legend([cfg.type ' ' cfg.DV], 'interpreter', 'none', 'autoupdate', 'off')
    %%
end % nGaits



%%
cd(cfg.datadir); cd ../../
%
cd(['Figures' filesep  cfg.type ' Fourier Fits'])
%%
print(['PFX overlayed ' cfg.type ' onset ' cfg.DV ' binned individual fourierfits'],'-dpng');
end % end function.