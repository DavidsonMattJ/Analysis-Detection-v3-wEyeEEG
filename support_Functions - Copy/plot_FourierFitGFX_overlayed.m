function plot_FourierFitGFX_overlayed(cfg, dataIN)




 figure(10); 
fntsize= 20;

GFX_headY = cfg.HeadData;

gaitfield = {'gc', 'doubgc'};
gaitprint = {'gait', 'stride'};
binfield = {'','_binned'};
nGaits_toPlot= cfg.nGaits_toPlot;
usebin=cfg.usebin;

barCols={'r', 'b', 'm'}'
pidx=cfg.pidx2;
mdiff = round(mean(diff(pidx)./2));
            xvec = pidx(1:end-1) + mdiff;
    
iLR=4; % doubgc?

    % plot each best fit.
    %collate data:
%     barCols= {'r',  'b', [.7 0 .7], [.2 .2 .2]};  
    p1=[];
    
    %which field of datastructure to plot?
    if strcmp(cfg.DV, 'RT')
        usefield = [gaitfield{cfg.nGaits_toPlot} binfield{cfg.usebin+1} '_rts'];        
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
    
    % extract data:
     
    %collate data:
    for isub= 1:size(dataIN,1)
        
        ppantData(isub,:)= dataIN(isub,iLR).(usefield);
         if nGaits_toPlot==2
            plotHead(isub,:) = GFX_headY(isub).([gaitfield{nGaits_toPlot} '_LRL']);
        else
            plotHead(isub,:) = GFX_headY(isub).(gaitfield{nGaits_toPlot});
        end
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
    
    
    
    % plot best fit (no bars though).
    gM = squeeze(mean(ppantData));
    
    
    [f,gof]= fit(xvec',gM',  'fourier1');
%     % % % % % % % %
%     subplot(2,1,2); 
    
    hold on;   
    h=plot(f, xvec, gM);%,
    h(2).LineWidth = 2;
    h(2).Color = barCols{cfg.iDV};
    h(1).Visible = 'off'; %remove data points.
    plotagain= h(2).YData;
    h(2).Visible = 'off';
    
    %replot but rescale, to show overlay:
    %rezero
    plotagain = plotagain - mean(plotagain);
    %rescale
    plotagain = plotagain ./ max(abs(plotagain));
%     if testtype==1
% %         plotagain= plotagain .*-1;
%     end
    p1(cfg.iDV) =plot(h(2).XData,  plotagain, 'color', barCols{cfg.iDV}, 'linew', 5, 'linest', '-');
    legend off
    set(gca,'fontsize', fntsize, 'xtick', [1, midp, xvec(end)], 'XTickLabels', {'0', '50', '100%'}, 'YTickLabels', {'-1', '-0.5', '0', '0.5', '1', ''})   
    xlabel([ 'position relative to ' gaitprint{nGaits_toPlot} '-cycle (%)']);%
    ylabel({['[a.u.]'];['(decrease  -  increase)']})
    
    ylim([-1 1.6])
    
    %% 

%head Y(doub)
if testtype==2
yyaxis right;
  pH= nanmean(plotHead,1); hold on;
% %   stEH = CousineauSEM(plotHead);
%     sh= shadedErrorBar(1:size(plotHead,2), pH, stEH,'k',1) ; 
%     sh.mainLine.LineStyle='-';
%     sh.mainLine.LineWidth = 4;
    ph = plot(1:size(plotHead,2), pH, 'color', 'k', 'linewidth', 4);
    legp = ph;
    set(gca,'ytick', [], 'fontsize',fntsize);
    
    midp=xvec(ceil(length(xvec)/2));
    lg=legend(legp, 'Head height [a.u.]', 'location', 'NorthEastOutside', 'fontsize', fntsize);
 set(gca,'fontsize', fntsize, 'xtick', [1, midp, xvec(end)], 'XTickLabels', {'0', '50', '100%'})   
    xlabel([ 'position relative to ' gaitprint{nGaits_toPlot} '-cycle (%)']);%
    ylim([-2 .6]);
    yyaxis left;
end

%%
legend([ph, p1(2), p1(1) p1(3)], ...
    {'head height','accuracy', 'reaction time', 'response onset'}, 'Location','NorthEastOutside');
end % plot summary (GFX fits overlayed)


