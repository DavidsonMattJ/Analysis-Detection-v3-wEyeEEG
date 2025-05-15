function  plot_PFfits_comparison_MSVer1_altFits(dataIN, cfg)
% Here we will plot a pretty MS version for FIgure 2 (first data figure).
% displaying condition (Stationary vs Walking) differences in performance.
%%

usecolsStWlk = {[.2 .2 .2], [.7 .7 0]}; % R Gr
usecolsLR= {[.7 0, 0], [0 .7 0], [.7, 0, .7]}; % R Gr Prple.
figure(1); clf; set(gcf, 'color', 'w', 'units', 'normalized', 'position', [0.05 0.05 .67 .65]);
nsubs = length(cfg.subjIDs);


fsize =15;
useLogscale=1;
if strcmp(cfg.plotlevel, 'GFX')
    clf;

    psubj= ['Group N=' num2str(nsubs)];
    % both this and the next use the same figure function:

    usegaitfields = {'gc_', 'doubgc_'};
    usetypes = {'_allstatnry', '_allwlking'};

   
    %extract data across ppants (for-loop once).
    acc_StWlk = zeros(nsubs,2); %mean acc
    rt_StWlk  =zeros(nsubs,2,7); % rt per contrast level.
    fitsObserved_StWlk =[];
    fitsObserved_Wlk_LR =[];
    
    [GFX_params_StWlk,GFX_params_Wlk_LR]=deal(zeros(nsubs,2,2)); % 2 parmas(threshold and slopes).
    GFX_contrastLevels = zeros(nsubs,2,7);
    % which plot type? only 1 gait supported for now.
    nGaits_toPlot=1;

    useg= usegaitfields{nGaits_toPlot};
    GFX_stimlists=[]; % for plotting the horizontal errorbars later.

    for ippant = 1:nsubs

        usetypes = {'_allstatnry', '_allwlking'};
        for iStWlk=1:2

%              for ipart=1:3
%                 %plot main and add extensions.
%                 pX=dataIN(ippant,3).([pfparts{ipart} 'X' usetypes{iStWlk}]);
%                 pY=dataIN(ippant,3).([pfparts{ipart} 'Y' usetypes{iStWlk}]);
% 
%                 plot(pX,pY, 'color', usecolsStWlk{iStWlk}, 'linestyle', uselns{ipart}, 'linew',2);
%                 hold on
%                 end
%                 % add markers:
%                 Stimlist=dataIN(ippant,3).(['data' usetypes{iStWlk}])(:,1);
%                 accper =  dataIN(ippant,3).(['data' usetypes{iStWlk}])(:,2)./dataIN(ippant,3).(['data' usetypes{iStWlk}])(:,3) ;
%                 lg(iStWlk)= plot(Stimlist, accper, '.', 'color', usecolsStWlk{iStWlk},'markersize',40');


StimList=dataIN(ippant,3).(['data' usetypes{iStWlk}])(:,1);
mAcc =  dataIN(ippant,3).(['data' usetypes{iStWlk}])(:,2)./dataIN(ippant,3).(['data' usetypes{iStWlk}])(:,3) ;
RTs=dataIN(ippant,3).(['data' usetypes{iStWlk}])(:,4);

% thresh, slope.
dt=dataIN(ippant,3).(['fitresult' usetypes{iStWlk}]);

% use the fullfit (incl extensions).
fullX = [dataIN(ippant,3).(['lowX' usetypes{iStWlk}]), dataIN(ippant,3).(['mainX' usetypes{iStWlk}]),...
dataIN(ippant,3).(['highX' usetypes{iStWlk}])];
fullY = [dataIN(ippant,3).(['lowY' usetypes{iStWlk}]), dataIN(ippant,3).(['mainY' usetypes{iStWlk}]),...
dataIN(ippant,3).(['highY' usetypes{iStWlk}])];

%>> store
acc_StWlk(ippant,iStWlk) = mean(mAcc);
acc_perContrast(ippant, iStWlk,:)= mAcc;

rt_StWlk(ippant, iStWlk,:) = RTs;
GFX_contrastLevels(ippant,iStWlk,:)=StimList;
GFX_params_StWlk(ippant, iStWlk,:) = [dt(1), dt(5)]; % els are Thresh, Width, Guess, Lapse , slope

fitsObserved_StWlk(ippant,iStWlk,1,:)= fullX;
fitsObserved_StWlk(ippant,iStWlk,2,:)= fullY;

        end
        %  now extract observed and fit per L/R ft.
     


        for iLR=1:2 % left / right step leading.

            StimList=dataIN(ippant,iLR).(['data_LRwlking'])(:,1);
            mAcc =  dataIN(ippant,iLR).(['data_LRwlking'])(:,2)./dataIN(ippant,iLR).(['data_LRwlking'])(:,3) ;
            RTs=dataIN(ippant,iLR).(['data_LRwlking'])(:,4);

            % thresh, slope.
            dt=dataIN(ippant,iLR).(['fitresult_LRwlking']);

            % use the fullfit (incl extensions).
            fullX = [dataIN(ippant,iLR).(['lowX_LRwlking']), dataIN(ippant,iLR).(['mainX_LRwlking']),...
                dataIN(ippant,iLR).(['highX_LRwlking'])];
            fullY = [dataIN(ippant,iLR).(['lowY_LRwlking']), dataIN(ippant,iLR).(['mainY_LRwlking']),...
                dataIN(ippant,iLR).(['highY_LRwlking'])];

            % >> store:

            fitsObserved_Wlk_LR(ippant,iStWlk,1,:)= fullX;
            fitsObserved_Wlk_LR(ippant,iStWlk,2,:)= fullY;
            GFX_params_Wlk_LR(ippant,iLR,:) = [dt(1), dt(5)]; % els are Thresh, Width, Guess, Lapse , slope
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
    ylabs={'Hit Rate', 'Reaction time [s]'};
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

        if accp>.05
            plotp = '\itns';
        else
            plotp= ['\itp \rm=' sprintf('%.3f', accp)];
        end
        text(mean(placeAt),yminmax(2)+.12*diff(yminmax) ,plotp, 'fontsize', fsize, 'HorizontalAlignment','center',...
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
        %revert to full search.
        ProportionCorrectObserved = squeeze(nanmean(acc_perContrast(:, iStWlk,:),1));

        hold on

        lg(iStWlk)=plot(StimList,ProportionCorrectObserved,['.'],'color', usecolsStWlk{iStWlk},'markersize',20);
        errObs = CousineauSEM(squeeze(acc_perContrast(:,iStWlk,:)));
        errorbar(StimList, ProportionCorrectObserved, errObs, 'k', 'linestyle', 'none');
        
        % add a group fit (mean of individual fits):
                av_IndFit_X = squeeze(mean(fitsObserved_StWlk(:,iStWlk,1,:),1));
        av_IndFit = squeeze(mean(fitsObserved_StWlk(:,iStWlk,2,:),1));

        err_IndFit = CousineauSEM(squeeze(fitsObserved_StWlk(:,iStWlk,2,:)));

%         plot(av_IndFit_X, av_IndFit, '-', 'color',usecolsStWlk{iStWlk} , 'linew', 2);
shadedErrorBar(av_IndFit_X, av_IndFit, err_IndFit, { 'linestyle', '-', 'color',usecolsStWlk{iStWlk} , 'linew', 2},1);

        if useLogscale==1 && iStWlk==1
            xlabel('Target contrast (dB)');
            xts= get(gca, 'XTickLabel');
            newl=[];
            for il= 1:length(xts);
                newl(il)= str2num(xts{il});
            end
            set(gca, 'XTickLabel', 20*log10(newl));
        elseif useLogscale==0
            xlabel('Target contrast');
        end
        ylabel(ylabs{1});        hold on
        %             title(psubj);
        xtickStim(iStWlk,:) = StimList;



        StimList_Group = StimList;
        StimList_Group_err = CousineauSEM(squeeze(GFX_contrastLevels(:,iStWlk,:)));


        qpos= .25+iStWlk*-.075; % where to place the qmarkers and errorbars?

        errorbar(StimList_Group, ProportionCorrectObserved, StimList_Group_err, ...
             'k', 'horizontal','linestyle','none', 'LineWidth',1);
        % add markers to clarify
        %
        %
        hghtMax = [.375 .5 .625 .75 .825 .855 .875]; % pos on PFit
        hghtMax= repmat(qpos+.025, [1, 7]);
%         Qn={'(-3)', '(-2)', '(-1)', 'QM', '(+1)', '(+2)', '(+3)'}
qbars={'(-3)', '(-2)', '(-1)', '\itM', '(+1)', '(+2)', '(+3)'};
%         Qn={'', '', '', 'QM', '', '', ''}

        offsetText= [.075, -.075];
        for iH = 1:length(StimList_Group)
            %draw vertical line to data.
            vertLims = [qpos hghtMax(iH)];
%             plot([StimList_Group(iH) StimList_Group(iH)], vertLims, ':', 'color', usecolsStWlk{iStWlk}, 'linew', 1)           
%             text(StimList_Group(iH), vertLims(1)-.03, Qn{iH},...
%                 'HorizontalAlignment', 'center','color', usecolsStWlk{iStWlk},...
%                 'fontsize', fsize-5, 'fontweight', 'bold');
%             
    indxAt = dsearchn(av_IndFit_X, StimList_Group(iH));

            text(StimList_Group(iH), av_IndFit(indxAt)+offsetText(iStWlk), qbars{iH},...
                'HorizontalAlignment', 'center','color', usecolsStWlk{iStWlk},...
                'fontsize', fsize-5);
        end


        % add threshold?
        ThreshG = dsearchn(av_IndFit, .75);

        plot([av_IndFit_X(ThreshG),av_IndFit_X(ThreshG)], [0 av_IndFit(ThreshG)], '-','color', usecolsStWlk{iStWlk});



    end

    %%
%     xlim([-33 -22])
xlim([0.415 .48])
    ylim([.1 1])
    shg

    %tidyax
    box off
    set(gca, 'fontsize',fsize)%, 'YTickLabel',[]);
    xlab= linspace(min(xtickStim(:)), max(xtickStim(:)),10);
    xlab= sort(xtickStim(:));
    xlab= xlab(1:2:length(xlab))
    if useLogscale==1

        set(gca, 'Xtick', xlab,  'xticklabel', round(20*log10(xlab),2));

    else
    set(gca, 'Xtick', xlab,  'xticklabel', round(xlab,2));
    end
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
qbars={'(-3)', '(-2)', '(-1)', '\itM', '(+1)', '(+2)', '(+3)'};
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    text(x,[.05,.05], qbars(i), 'HorizontalAlignment', 'center', 'fontsize', 8);
end



    ylabel('Reaction time [s]');
    xlabel('')
    ylim([0 .5]);
%     text(1.5,.4375,['\itp \rm < .001'], 'fontsize', fsize, 'HorizontalAlignment','center');
%     plot([1,2], [.425 .425], 'k-', 'linew', 2)
shg
box off

   
    %%
    cd([cfg.figdir  filesep 'Accuracy summary'])
    print('-dpng', [psubj ' accuracy and rt summary MS ver1 altFITS']);

% can also push to JASP for rmANOVAs:
% 
% GFX_rt_StWlk=table();
% wlks={'stat', 'wlk'};
% for iwlk= 1:2
%     for iQ=1:7
% GFX_rt_StWlk.([wlks{iwlk} '_' num2str(iQ)]) = rt_StWlk(:,iwlk,iQ);
% 
%     end
% end
% cd ../
% %write for jasp
% writetable(GFX_rt_StWlk, 'GFX_rt_Stwlk.csv');


end% GFX


end