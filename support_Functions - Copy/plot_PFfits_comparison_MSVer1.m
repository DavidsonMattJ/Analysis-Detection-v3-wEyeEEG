function  plot_PFfits_comparison_MSVer1(dataIN, cfg)
% Here we will plot a pretty MS version for FIgure 2 (first data figure).
% displaying condition (Stationary vs Walking) differences in performance.
%%

GFX_headY = cfg.HeadData;
usecolsStWlk = {[.2 .2 .2], [.7 .7 0]}; % R Gr
usecolsLR= {[.7 0, 0], [0 .7 0], [.7, 0, .7]}; % R Gr Prple.
figure(1); clf;
set(gcf, 'color', 'w', 'units', 'normalized', 'position', [0.05 0.05 .67 .75]);
nsubs = length(cfg.subjIDs);


fsize =15;
useLogscale=1;

if useLogscale==1
    PF = @PAL_Logistic; % if on log scale. (stim intensities are linear for now)
else
    % PF= @PAL_Quick; %if stim are linearly transformed
    PF= @PAL_Weibull; %if stim are linearly transformed
end

if strcmp(cfg.plotlevel, 'GFX')
    clf;

    psubj= ['Group N=' num2str(nsubs)];
    % both this and the next use the same figure function:

    usegaitfields = {'gc_', 'doubgc_'};
    usetypes = {'_allstatnry', '_allwlking'};

    %% Palamedes parameter grid, used across all plots:
    % %Parameter grid defining parameter space through which to perform a
    % %brute-force search for values to be used as initial guesses in iterative
    % %parameter search.


    searchGrid.beta = logspace(0,3,101); % slope grid
    %         searchGrid.gamma = 0:.01:.15; % guess rate
    %         searchGrid.lambda= 0:.01:.15; % lapse rate
    searchGrid.gamma = 0.1;  %guess rate, scalar here (since fixed) but may be vector
    searchGrid.lambda = 0.02;  %lapse rate

    %Threshold and Slope are free parameters, guess and lapse rate are fixed
    paramsFree = [1 1 0 0];  %1: free parameter, 0: fixed parameter

    %extract data across ppants (for-loop once).
    acc_StWlk = zeros(nsubs,2); %mean acc
    rt_StWlk  =zeros(nsubs,2,7); % rt per contrast level.
    fitsObserved_StWlk = zeros(nsubs,2,7); %acc per contrast level.

    fitsObserved_Wlk_LR = zeros(nsubs,2,7); %fits split by leading foot.
    [GFX_params_StWlk,GFX_params_Wlk_LR]=deal(zeros(nsubs,2,2)); % 2 parmas(threshold and slopes).
    GFX_contrastLevels = zeros(nsubs,2,7);
    % which plot type? only 1 gait supported for now.
    nGaits_toPlot=1;

    useg= usegaitfields{nGaits_toPlot};
    GFX_stimlists=[]; % for plotting the horizontal errorbars later.
    for ippant = 1:nsubs

        usetypes = {'_allstatnry', '_allwlking'};
        for iStWlk=1:2

            StimList= dataIN(ippant,1).(['gc_ContrVals' usetypes{iStWlk}]);

            if useLogscale==1
                StimList= StimList- .4; % relative to background
                % note contrast was on range of .4:1;
                StimList = 20*log10(StimList);
            end
            StimLevelsFineGrain=linspace(min(StimList),max(StimList),200);

            searchGrid.alpha= (StimLevelsFineGrain);

            GFX_contrastLevels(ippant,iStWlk,:)=( StimList);


            NumPer = dataIN(ippant,1).([useg 'NumPerContr' usetypes{iStWlk}]);
            TotalPer = dataIN(ippant,1).([useg 'TotalPerContr' usetypes{iStWlk}]);
            RTs = dataIN(ippant,1).([useg 'RTperContr' usetypes{iStWlk}]);
            % take mean accuracy over all contrast vals.
            mAcc = sum(NumPer) / sum(TotalPer);

            acc_StWlk(ippant,iStWlk) = mAcc;

            rt_StWlk(ippant, iStWlk,:) = RTs;

            % define guess rate relativey to minimum participant accuracy.
            searchGrid.gamma = min(NumPer./TotalPer)-.05;

            [paramsValues, LL, scenario,output] = PAL_PFML_Fit(StimList,NumPer, ...
                TotalPer,searchGrid,paramsFree,PF, 'lapseLimits', [0 1],'guessLimits',...
                [0 1]);

            if paramsValues(1) < min(StimList)
                disp('check')
            end
            fitsObserved_StWlk(ippant,iStWlk,:) = NumPer./ TotalPer;
            %                 fitsModel_StWlk(ippant).d(iStWlk,:) = PF(paramsValues,StimLevelsFineGrain);

            GFX_params_StWlk(ippant, iStWlk,:) = paramsValues(1:2);

            % participant summary
            disp(['Threshold est, ppant ' num2str(ippant) ': '])

            GFX_stimlists(ippant,iStWlk,:) = StimList;
        end


        % extract observed and fit per L/R ft.
        usetypes = {'_LRwlking'};



        for iLR=1:2 % left / right step leading.

            NumPer = dataIN(ippant,iLR).([useg 'NumPerContr_LRwlking']);
            TotalPer = dataIN(ippant,iLR).([useg 'TotalPerContr_LRwlking']);
            NumPer(NumPer==0)= .001;

            StimList = dataIN(ippant,iLR).([useg 'ContrVals_LRwlking']);
            if useLogscale==1
                StimList= 20*log10(StimList);
            end
            contrVals_grain = linspace(min(StimList),max(StimList),200);
            searchGrid.alpha = contrVals_grain;

            [paramsValues LL exitflag] = PAL_PFML_Fit(StimList,NumPer, ...
                TotalPer,searchGrid,paramsFree,PF);


            fitsObserved_Wlk_LR(ippant,iLR,:) = NumPer./TotalPer;
            %                 fitsModel_Wlk_LR(ippant,iLR,:) =  PF(paramsValues,StimLevelsFineGrain);

            GFX_params_Wlk_LR(ippant,iLR,:) = paramsValues(1:2);
        end % LR


    end % ippant

    %% show contrast levels per ppant.
    %         figure(9);
    %         plot(1:size(GFX_contrastLevels,1), GFX_contrastLevels', 'o', 'color','b')
    %         ylabel('contrast'); xlabel('ppant');



    %% % continue with plots:
    %% First grand effects by condition (stationary vs walking) %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    legp=[]; % for legend

    useD= {acc_StWlk, squeeze(mean(rt_StWlk,3))};
    pspots = [1,6] ; % subplot position.
    ylabs={'Accuracy', 'Reaction time [s]'};
    storeYlims=[];
    for id=1:2
        tmpD= useD{id};
        mBar = squeeze(nanmean(tmpD,1));
        errBar = CousineauSEM(tmpD);
        %%
        subplot(2,5,pspots(id));
        cla; % thin plot.
        placeAt = [.7 2.3];
        %         bh1= bar([mBar(1), nan], 'FaceAlpha', .5, 'FaceColor', usecolsStWlk{1} );            hold on
        %         bh2= bar([ nan, mBar(2)], 'FaceAlpha', .5, 'FaceColor', usecolsStWlk{2});
        ph1 = plot(placeAt(1), mBar(1), 'o','Color', usecolsStWlk{1}, 'MarkerSize', 20, 'LineWidth',2); hold on;
        ph1 = plot(placeAt(2), mBar(2), 'o', 'Color', usecolsStWlk{2}, 'MarkerSize',20, 'LineWidth',2);
        shg
        xlim([.5 2.5])
        %
        errorbar(placeAt, mBar, errBar, 'k','linestyle', 'none', 'linew', 2);
        %include individual dpoints:
        scD= ones(1,nsubs);
        %add jitter to display ind subj points:
        jit = (rand(1,100) - .5)./8;
        plot(1+jit(1:nsubs), tmpD(:,1),'.', 'color',  [usecolsStWlk{1}, .9], 'markersize', 20)
        plot(2+jit(1:nsubs), tmpD(:,2),'.', 'color',  [usecolsStWlk{2}, .9], 'markersize', 20)
        % connect per subj.
        for ippant= 1:nsubs
            plot([1+jit(ippant) 2+jit(ippant)], [tmpD(ippant,1), tmpD(ippant,2)], 'color', [.8 .8 .8], 'linew', 1)
        end

        %tidy axes
        ylabel(ylabs{id});
        %         title(psubj);
        set(gca, 'Xtick', placeAt,'xticklabels', {'standing', 'walking'}, 'fontsize', fsize);
        % add sig:
        %add text with relv stats:
        [~,accp, ~, stats] = ttest(tmpD(:,1),  tmpD(:,2));

        d=computeCohen_d(tmpD(:,1),  tmpD(:,2),'paired');
        xlim([0 3]);
        %define y lim based on width of points.
        yminmax=[ min(tmpD(:)), max(tmpD(:))];

        % ylim from min:max +25%
        ylim([yminmax(1)- .25*(diff(yminmax)), yminmax(2)+ .3*(diff(yminmax))])

        storeYlims(id,:) = get(gca,'Ylim');
        text(mean(placeAt),yminmax(2)+.12*diff(yminmax) ,['\itp \rm=' sprintf('%.3f', accp)], 'fontsize', fsize, 'HorizontalAlignment','center',...
            'VerticalAlignment','bottom');

        plot(placeAt, [yminmax(2)+.1*diff(yminmax) yminmax(2)+.1*diff(yminmax)], 'k-', 'linew', 2)
        box off
    end

    %% Now PF FITS %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %now for PF fits.
    lg=[];
    %for group level plots, change x axis to 1:7
    xtickStim=[];
    for iStWlk=1:2

        subplot(2,5,[2:5]);
        % Plot the group fits, but do stats at individual level
        % (compare paramValues)
        %
        %             ProportionCorrectModel =squeeze(nanmean(fitsModel_StWlk(:,iStWlk,:),1));
        %             errModel = CousineauSEM(squeeze(fitsModel_StWlk(:,iStWlk,:)));
        %             shadedErrorBar(StimLevelsFineGrain, ProportionCorrectModel, errModel, {'color', usecolsStWlk{iStWlk}});
        %
        StimList= squeeze(mean(GFX_contrastLevels(:,iStWlk,:),1)); % already in linear or log
        contrVals_grain = (linspace(min(StimList),max(StimList),200));
        StimList = (StimList);
        searchGrid.alpha = (contrVals_grain);
        %revert to full search.
        searchGrid.gamma = 0.05;  %guess rate, scalar here (since fixed) but may be vector
        searchGrid.lambda = 0.02;  %lapse rate
        ProportionCorrectObserved = squeeze(nanmean(fitsObserved_StWlk(:, iStWlk,:),1));
        hold on

        lg(iStWlk)=plot(StimList,ProportionCorrectObserved,['.'],'color', usecolsStWlk{iStWlk},'markersize',20);
        errObs = CousineauSEM(squeeze(fitsObserved_StWlk(:,iStWlk,:)));
        errorbar(StimList, ProportionCorrectObserved, errObs, 'k', 'linestyle', 'none');
        %%
        %%add group fit:
        %             [paramsValues LL exitflag] = PAL_PFML_Fit(StimList,round(ProportionCorrectObserved.*100)', ...
        %                 repmat(100,1, length(ProportionCorrectObserved)),searchGrid,paramsFree,PF);
        paramsValues = PAL_PFML_Fit(round(StimList,3)',ProportionCorrectObserved', ...
            ones(1, length(ProportionCorrectObserved)),searchGrid,paramsFree,PF);

        groupFit = PF(paramsValues,contrVals_grain);

        plot(contrVals_grain, groupFit, '-', 'color',usecolsStWlk{iStWlk} , 'linew', 2);

        if useLogscale==1
            xlabel('Target contrast (dB)');
        else
            xlabel('Target contrast');
        end
        ylabel('Accuracy');
        hold on
        %             title(psubj);
        xtickStim(iStWlk,:) = StimList;



        StimList_Group = squeeze(nanmean(GFX_stimlists(:, iStWlk,:),1));
        StimList_Group_err = CousineauSEM(squeeze(GFX_stimlists(:, iStWlk,:)));


        qpos= .35+iStWlk*-.075; % where to place the qmarkers and errorbars?


        % plot the stim lists (horizontal errorbars)
%         errorbar(StimList_Group, repmat(qpos, 1, length(StimList_Group)),...
%             StimList_Group_err, 'k', 'horizontal','linestyle','none', 'LineWidth',1)
%         
        
        errorbar(StimList_Group, ProportionCorrectObserved, StimList_Group_err, ...
             'k', 'horizontal','linestyle','none', 'LineWidth',1);
        % add markers to clarify
        %
        %
        hghtMax = [.375 .5 .625 .75 .825 .855 .875]; % pos on PFit
%         Qn={'(-3)', '(-2)', '(-1)', 'QM', '(+1)', '(+2)', '(+3)'}
        Qn={'Q(-3)', 'Q(-2)', 'Q(-1)', 'QM', 'Q(+1)', 'Q(+2)', 'Q(+3)'}
%         Qn={'', '', '', 'QM', '', '', ''}

        for iH = 1:length(StimList_Group)
            %draw vertical line to data.
            vertLims = [qpos hghtMax(iH)];
            plot([StimList_Group(iH) StimList_Group(iH)], vertLims, ':', 'color', usecolsStWlk{iStWlk}, 'linew', 1)
            text(StimList_Group(iH), vertLims(1)-.03, Qn{iH},...
                'HorizontalAlignment', 'center','color', usecolsStWlk{iStWlk},...
                'fontsize', fsize-5, 'fontweight', 'bold');
        end

    end
    xlim([-33 -22])
    ylim([.1 1])
    shg

    %tidyax
    box on
    set(gca, 'fontsize',fsize)%, 'YTickLabel',[]);
    xlab= linspace(min(xtickStim(:)), max(xtickStim(:)),10);
    set(gca, 'Xtick', xlab,  'xticklabel', round(xlab,2));

    legend([lg], {'standing', 'walking'}, 'location', 'NorthWest', 'fontsize',fsize, 'autoupdate', 'off' )
    %add text with relv stats:
    [~,threshsig, ~, thrstats] = ttest(squeeze(GFX_params_StWlk(:,1,1)),  squeeze(GFX_params_StWlk(:,2,1)));
    [~,slopesig, ~, slpstats] = ttest(squeeze(GFX_params_StWlk(:,1,2)),  squeeze(GFX_params_StWlk(:,2,2)));

    threshd = computeCohen_d(squeeze(GFX_params_StWlk(:,1,1)),  squeeze(GFX_params_StWlk(:,2,1)));

    sloped = computeCohen_d(squeeze(GFX_params_StWlk(:,1,2)),  squeeze(GFX_params_StWlk(:,2,2)));

    textmsg = {['thresh. \itt\rm(' num2str(thrstats.df) ')=' sprintf('%.2f',thrstats.tstat) ', \itp \rm= ' sprintf('%.3f', threshsig) ];...
        ['slope \itt\rm(' num2str(slpstats.df) ')=' sprintf('%.2f',slpstats.tstat) ', \itp \rm= ' sprintf('%.3f', slopesig) ]};

    %         text(max(StimList),.4,textmsg, 'fontsize', 12, 'HorizontalAlignment', 'right');
    disp(textmsg)


    %% now RTs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    subplot(2,5,7:10); cla;
    mBar =  squeeze(nanmean(rt_StWlk,1));
    st1 = CousineauSEM(squeeze(rt_StWlk(:,1,:)));
    st2 = CousineauSEM(squeeze(rt_StWlk(:,2,:)));
    errBar = [st1;st2];
    % 
        bh1= bar([mBar(1,:); nan(1,7)], 'FaceAlpha', .5, 'FaceColor', usecolsStWlk{1});            hold on
        bh2= bar([ nan(1,7); mBar(2,:)], 'FaceAlpha', .5, 'FaceColor', usecolsStWlk{2});
%         legend([bh1(1), bh2(2)], {'standing', 'walking'}, 'autoupdate', 'off')
        %
        errorbar_groupedfit( mBar, errBar);
%       set(gca,'fontsize',fsize,'xtick', [1:7], 'xticklabel', {'Q(-3)', 'Q(-2)', 'Q(-1)', 'QM', 'Q(+1)', 'Q(+2)', 'Q(+3)'});
    set(gca,'fontsize',fsize,'xtick', [1:2], 'xticklabel', {'standing', 'walking'});

    % add text to bars:
    [ngroups, nbars]= size(mBar);
    groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
qbars={'Q(-3)', 'Q(-2)', 'Q(-1)', 'QM', 'Q(+1)', 'Q(+2)', 'Q(+3)'};
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    text(x,[.05,.05], qbars(i), 'HorizontalAlignment', 'center', 'fontsize', fsize-5);
end



    ylabel('Reaction time [s]');
    xlabel('')
    ylim([0 .5]);
%     text(1.5,.4375,['\itp \rm < .001'], 'fontsize', fsize, 'HorizontalAlignment','center');
%     plot([1,2], [.425 .425], 'k-', 'linew', 2)
shg


   
    %%
    cd([cfg.figdir  filesep 'Accuracy summary'])
    print('-dpng', [psubj ' accuracy and rt summary MS ver1']);

% can also push to JASP for rmANOVAs:

GFX_rt_StWlk=table();
wlks={'stat', 'wlk'};
for iwlk= 1:2
    for iQ=1:7
GFX_rt_StWlk.([wlks{iwlk} '_' num2str(iQ)]) = rt_StWlk(:,iwlk,iQ);

    end
end
cd ../
%write for jasp
writetable(GFX_rt_StWlk, 'GFX_rt_Stwlk.csv');


end% GFX


end