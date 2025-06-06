function plot_PFfits_gaitcycle_GFX_MSver1(dataIN, cfg)

% helper function to plot the psychometric fits at the Group level, showing
% difference over the gait cycle.
%based on plot_PFfits_gaitcycle.m

GFX_headY = cfg.HeadData;
figure(1); clf; set(gcf, 'color', 'w', 'units', 'normalized', 'position', [0 0 .9 .9]);
nsubs = length(cfg.subjIDs);

%set colours for gait quintiles.
% qntlCols = cbrewer('div', 'Spectral', 100);
% qntlCols = cbrewer('seq', 'Reds', 100);
qntl1 = cbrewer('seq','Reds',50);
qntl2 = flipud(cbrewer('seq','Blues',50));
qntlCols=[qntl1;qntl2];
qntlCols(qntlCols>1)=1;

% load a circular colormap. romO bamO brocO corkO vikO
load('vikO.mat');
qntlCols = vikO;
%resampe and adjust.
qntlCols= imresize(qntlCols, [100,3]);
useLogscale=1;

if useLogscale==1
    PF = @PAL_Logistic; % if on log scale. (stim intensities are linear for now)
else
    % PF= @PAL_Quick; %if stim are linearly transformed
    PF= @PAL_Weibull; %if stim are linearly transformed
end


plotGroupLevel_Pfit=2; % 1 for Psych fit assessed on group level, 2 for average of ppant fits.


if strcmp(cfg.plotlevel, 'GFX')
    clf;
    
    psubj= 'GFX';
    % both this and the next use the same figure function:
    
    usegaitfields = {'gc', 'doubgc'};
    usetypes = {'_qntlwlking'};
    %% Palamedes parameter grid, used across all plots:
    % %Parameter grid defining parameter space through which to perform a
    % %brute-force search for values to be used as initial guesses in iterative
    % %parameter search.
    
    %contrast values across exp (perppant):
   
    searchGrid.beta = logspace(0,3,101);
    searchGrid.gamma = 0.1;  %scalar here (since fixed) but may be vector
    searchGrid.lambda = 0.02;  %ditto
    
    %Threshold and Slope are free parameters, guess and lapse rate are fixed
    paramsFree = [1 1 0 0];  %1: free parameter, 0: fixed parameter
   
    
    %extract data across ppants (for-loop once).
    GFX_headdata= zeros(nsubs, 100);
    ppantfitsObserved_WlkQuantile = zeros(nsubs, 2,3, 7); % contrast levels.
    ppantfitsModel_WlkQuantile = zeros(nsubs,2, 3,200);
    ppantfitsStimRange_WlkQuantile= zeros(nsubs,2, 3,200);
    GFX_params_WlkQuantile=zeros(nsubs,2,3,2); % subs, iLR, quantiles, slope/thresh.
    GFX_stimlists= zeros(nsubs,3,7); % subs, iLR, 7 levels.
    for nGaits_toPlot=1%:2 % 2 not yet supported
        
        useg= usegaitfields{nGaits_toPlot};
        
        qntlBounds = dataIN(1,1).([useg '_qntl_bounds']); % should be the same across ppants.
    
        for ippant=1:nsubs
            GFX_headdata(ippant,:) = GFX_headY(ippant).gc;
            %contrast values across exp (perppant):
          
            % note that other searchGrid params set above.

             for iLR=1:3 % 3rd dim is combined.
                 NumPerQ= dataIN(ippant, iLR).([useg '_NumPerContr_qntlwlking']);
                 TotalPerQ = dataIN(ippant,iLR).([useg '_TotalPerContr_qntlwlking']);
                 StimList= dataIN(ippant,iLR).gc_ContrVals_allwlking;
                 StimLevelsFineGrain=linspace(min(StimList),max(StimList),200);

                 if useLogscale==1
                     StimList= 10*log10(StimList);
                     StimLevelsFineGrain=10* log10(StimLevelsFineGrain);
                 end
            searchGrid.alpha = StimLevelsFineGrain;

                for iqnt= 1:size(NumPerQ,1)
                    %qntl specific data:
                    NumPer = NumPerQ(iqnt,:);
                    TotalPer = TotalPerQ(iqnt,:);
                    %% remove any zeros:
                    NumPer(NumPer==0)= .001;
                    [paramsValues LL scenario output] = PAL_PFML_Fit(StimList,NumPer, ...
                        TotalPer,searchGrid,paramsFree,PF);
                    
                    %store observed data and ppant fit:
                    ppantfitsObserved_WlkQuantile(ippant, iLR,iqnt,:)=NumPer./TotalPer;
                    
                    if any(isinf(paramsValues))
                        disp('COULD NOT FIT PFIT to this ppant');
                        paramsValues(1:2)=nan(1,2);
                    end
                    
                GFX_params_WlkQuantile(ippant,iLR,iqnt,:) = paramsValues(1:2);
                    
                %also store the evaluated ppant fit. (alternate plot).
                ppantFit = PF(paramsValues,StimLevelsFineGrain);
                ppantfitsModel_WlkQuantile(ippant,iLR,iqnt,:) = ppantFit;
                ppantfitsStimRange_WlkQuantile(ippant,iLR,iqnt,:) = StimLevelsFineGrain;
                end % iqnt
                GFX_stimlists(ippant,iLR,:) = StimList;
            end % iLR
        end % ippant.
        
        %% % continue with plots: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         ftnames= {'LR', 'RL'};
        ftnames = {'L foot -> R foot', 'R foot -> L foot', 'all'};
    
%%
        for iLR=3
%             subplot(3,3,1 + (3*(iLR-1)));
            % plot Head data
subplot(2,3,1)
            mP= nanmean(GFX_headdata,1);
%             errP = CousineauSEM(GFX_headdata);
            plot(1:size(GFX_headdata,2),mP, ['k-o'])
            yls = get(gca, 'ylim');            
            % place patch color in BG.
            for iq = 1:length(qntlBounds)-1
                tidx = qntlBounds(iq):qntlBounds(iq+1);
                hold on;
                thisColr= qntlCols(qntlBounds(iq),:);
                plot(tidx, mP(tidx), 'color',thisColr , 'linew', 3);
                %add patch ?
                
                xv = [tidx(1) tidx(1) tidx(end) tidx(end)];
                yv = [yls(1) yls(2) yls(2) yls(1)];
                ph = patch('XData', xv, 'YData', yv, 'FaceColor', thisColr, 'FaceAlpha', .6 , 'LineStyle','none');
                %             ph.FaceColor = thisColr;
            end
            ylim(yls);
            title([ psubj ' ' ftnames{iLR}]);
            set(gca, 'fontsize', 12, 'ytick', []);
            ylabel('head height')
            xlabel('% step completion')
           
            
            %% now plot the psychometric curves for each portion.
          StimList_Group = squeeze(nanmean(GFX_stimlists(:, iLR,:),1))';
          StimList_Group_err = CousineauSEM(squeeze(GFX_stimlists(:, iLR,:)));
          StimLevelsFineGrain_Group=linspace(min(StimList_Group),max(StimList_Group),200);

                      searchGrid.alpha = StimLevelsFineGrain_Group;
                      searchGrid.beta = logspace(0,3,101);
                      searchGrid.gamma = 0.1;  %scalar here (since fixed) but may be vector
                      searchGrid.lambda = 0.02;  %ditto
%%
adjustX=[-.03:.01:.03]
jitwth = .02;
% cla
legh=[];
            for iqnt= 1:size(NumPerQ,1)
            
                %per quantile, add a small amount of jitter.
%                 jitt = (rand(1)-.5)/10;
                subplot(2,3,[2,3])
                hold on;
%                 plotX = StimList_Group+adjustX(iqnt);
jit = (rand(1,7) - 0.5) * jitwth;

                    plotX= StimList_Group + jit;
                   thisColr= qntlCols(qntlBounds(iqnt),:);
               %plot mean observed values per contrast, with errorbarS:                
                ProportionCorrectObserved = squeeze(nanmean(ppantfitsObserved_WlkQuantile(:, iLR,iqnt,:),1));
               %plot data:
                plot(plotX,ProportionCorrectObserved,['o'],'color','k','markersize',10);  
                pl=plot(plotX,ProportionCorrectObserved,['o'],'color',thisColr,'markersize',8);  
                pl.MarkerFaceColor = thisColr;
%                 pl.Ed   
                lg(iLR)=pl;
                %plot error (Vertical error in proportion correct):
                  errObs = CousineauSEM(squeeze(ppantfitsObserved_WlkQuantile(:,iLR,iqnt,:)));                
                  errorbar(plotX, ProportionCorrectObserved, errObs, 'k', 'linestyle', 'none');
%add group fit:
                [paramsValues LL exitflag] = PAL_PFML_Fit(plotX,round(ProportionCorrectObserved.*100)', ...
                             repmat(100,1, length(ProportionCorrectObserved)),searchGrid,paramsFree,PF);
%                 plot GROUP assessed fit, or mean of ppant fits?
                if plotGroupLevel_Pfit==1
                groupFit = PF(paramsValues,StimLevelsFineGrain_Group); 
                plot(StimLevelsFineGrain_Group, groupFit, '-', 'color',thisColr , 'linew', 2);

                else
                groupFit = squeeze(nanmean(ppantfitsModel_WlkQuantile(:,iLR,iqnt,:),1));
                stE= CousineauSEM(squeeze(ppantfitsModel_WlkQuantile(:,iLR,iqnt,:)));
%                 sh=shadedErrorBar(StimLevelsFineGrain_Group, groupFit, stE, {'color', thisColr},1);
%                 sh.edge(1).Color = 'k';
%                 sh.edge(2).Color = 'k';
ph=plot(StimLevelsFineGrain_Group, groupFit, 'Color', thisColr, 'LineWidth', 2);
                legh(iqnt)=ph;
                end

    set(gca,'color', 'w')
            end % iqnt
%
            %plot error (horizontal error, spread of ppant bins)
                    errorbar(StimList_Group, repmat(.35, 1, length(StimList_Group)),...
                        StimList_Group_err, 'k', 'horizontal','linestyle','none', 'LineWidth',1)
% add markers to clarify
%
width = [-.05*jitwth, .05 *jitwth];
hghtMax = [.375 .5 .625 .75 .825 .855 .875];
for iH = 1:length(StimList_Group)
%draw vertical line to data.
vertLims = [.35 hghtMax(iH)];
plot([StimList_Group(iH) StimList_Group(iH)], vertLims, 'k:')

end
%
legend(legh, {'qntl1', 'qntl2', 'qntl3', 'qntl4', 'qntl5'}, 'location',' NorthWest');

%%
             %add text with relv stats: (need rmANOVA)
        % convert data to table:
        subjs = [1:nsubs]';
        [Fs,Ps]=deal(zeros(2,1));
        
        for iThrSl=1:2 % for slope and threshold
            %make columns of the repeated measures (each quantile).
            measures=  squeeze(GFX_params_WlkQuantile(:, iLR, :,iThrSl));       
             t=splitvars(table(subjs, measures));
             %to make it robust to toggling gait quantiles, extract the
             %autonames:
             %wilkinson notation for a rmmodel.
             rmmodel = ['measures_1-measures_' num2str(length(qntlBounds)-1) '~1'];
        
             WithinDesign = table([1:size(measures,2)]','VariableNames',{'quantiles'});

             mdlfit= fitrm(t, rmmodel, 'WithinDesign', WithinDesign);
             
             %get stats from table:
             rtable = ranova(mdlfit);
             Fs(iThrSl,1) = rtable.F(1);
             Ps(iThrSl,1) = rtable.pValue(1);
        end
             
        textmsg1 = ['\itF\rm(' num2str(rtable.DF(1)) ',' num2str(rtable.DF(2)) ')=' sprintf('%.2f',Fs(1)) ', \itp\rm=' sprintf('%.2f', Ps(1))];
        textmsg2 = ['\itF\rm(' num2str(rtable.DF(1)) ',' num2str(rtable.DF(2)) ')=' sprintf('%.2f',Fs(2)) ', \itp\rm=' sprintf('%.2f', Ps(2))];
        text(4,.2,{textmsg1; textmsg2}, 'fontsize', 10);
        
%             xlim([.5 7.5])
            xlabel('contrast')
           
            if iLR==1
                title('fits over the gait');
            end
%             set(gca, 'Xtick',StimList_Group, 'xticklabel', split(num2str(StimList-4))', 'color', [.8 .8 .8]);
            
        if useLogscale==1
       xlabel('contrast increment (dB)');       
        else
        xlabel('contrast increment');
        end
        ylabel('proportion correct');
        set(gca, 'fontsize', 12);
        
        
        %% add summary bar chart to clarify effect:
        % grouped bar chart:
        barDataThresh = squeeze(GFX_params_WlkQuantile(:, iLR, :,1));
        barDataSlope = squeeze(GFX_params_WlkQuantile(:, iLR, :,2));   
        
        barDataPlot = {barDataThresh, barDataSlope};
        ytitles = {'threshold values', 'slope parameter'};
        ttitles = {['\rm' textmsg1], ['\rm' textmsg2]};
        subplotspot= [5,6];
        ylimsare=[0,.6; 0,200];
        
        for ithrSlp = 1:2 % 5 and 6 or subplot pos.
          subplot(2,2,2+ithrSlp);
          cla
        barD = barDataPlot{ithrSlp};
           mBar = nanmean(barD,1);
       bh= bar(1:size(barD,2),mBar,'FaceColor','flat'); %flat allows us to update colours below
       bh.CData = qntlCols(qntlBounds(1:end-1),:);
       hold on;
       stE = CousineauSEM(barD);
       errorbar(1:size(barD,2),mBar, stE,'linestyle', 'none')
%        ylim(ylimsare(ithrSlp,:))
% adjust ylims to sdrange.
        sdrange = max(mBar) - min(mBar);
        ylim([min(mBar)-1.5*sdrange max(mBar)+1.5*sdrange])
       title(ttitles{ithrSlp})
       ylabel(ytitles{ithrSlp})
       xlabel('quantile');
       set(gca,'fontsize',12)
%        axis tight
        end % each param to plot.

%%

        end % iLR
       
    end % nGaits
        cd([cfg.figdir filesep  cfg.type ' onset Qntl PsychFits'])
        
        print([psubj ' ' cfg.type ' onset qntl PsychFits_MS'],'-dpng');
        
end% GFX
end