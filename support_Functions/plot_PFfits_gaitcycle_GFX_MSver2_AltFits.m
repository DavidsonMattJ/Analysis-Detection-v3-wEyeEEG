function plot_PFfits_gaitcycle_GFX_MSver2_AltFits(dataIN, cfg)

% helper function to plot the psychometric fits at the Group level, showing
% difference over the gait cycle.
%based on plot_PFfits_gaitcycle.m
%%
GFX_headY = cfg.HeadData;
figure(1); clf; 
set(gcf, 'color', 'w', 'units', 'normalized', 'position', [0 0 .7 .7]);
nsubs = length(cfg.subjIDs);

%set colours for gait quintiles.
% qntlCols = cbrewer('div', 'Spectral', 100);
% qntlCols = cbrewer('seq', 'Reds', 100);
qntl1 = cbrewer('seq','Reds',50);
qntl2 = flipud(cbrewer('seq','Blues',50));
qntlCols=[qntl1;qntl2];
qntlCols(qntlCols>1)=1;

% load a circular colormap. romO bamO brocO corkO vikO
load('vikO.mat'); %red blue
load('romaO.mat'); %like jet
load('bamO.mat'); %pink purple green
load('brocO.mat'); %blue brown yellow
load("corkO.mat"); %blue green grey
qntlCols = vikO;
%resampe and adjust.
qntlCols= circshift(imresize(qntlCols, [100,3]),25);

useLogscale=1;

debugPlots=cfg.debug;
FitatGroup=cfg.FitatGroup;
shareX= cfg.shareX; % share the contrast bins for quintiles.
cleanOnly=cfg.cleanOnly; % restrict to subset of pparticipants with no outliers in any quartiles?
applyOffset = cfg.applyOffset;

fitstoUse= cfg.FitsfromSharedX+1;
extractFits= {'', '_sharedX'}; % so if cfg is set to 1, we extract the shared fit for plots).
if strcmp(cfg.plotlevel, 'GFX')
    clf;
    
    psubj= 'GFX';
    % both this and the next use the same figure function:

   
    %extract data across ppants (for-loop once).
    GFX_headdata= zeros(nsubs, 100);
  
    GFX_fitsObserved_perQntl=[];
    GFX_stimlists= zeros(nsubs,5,7); % subs, qntls, 7 levels.
    GFX_params_qntl=zeros(nsubs,5,3); % thresh, width, slope
    GFX_Acc_perQntl= zeros(nsubs, 5, 7); % qntls and contrasts.
    qntlBounds= [1,20,40,60,80,100];
    
        for ippant=1:nsubs
            GFX_headdata(ippant,:) = GFX_headY(ippant).gc;
          
            for iqntl= 1:5

                extractfield = ['q' num2str(iqntl) extractFits{fitstoUse}];

                NumPerQ= dataIN(ippant, iqntl).(['data_' extractfield])(:,2);
                TotalPerQ= dataIN(ippant, iqntl).(['data_' extractfield])(:,3);
                StimList= dataIN(ippant, iqntl).(['data_' extractfield])(:,1);

                mAcc= NumPerQ./TotalPerQ;
               
                % thresh, slope.
                dt=dataIN(ippant,iqntl).(['fitresult_' extractfield]);
%                 dt.Fit
                % use the fullfit (incl extensions).
                fullX = [dataIN(ippant,iqntl).(['lowX_' extractfield]), dataIN(ippant,iqntl).(['mainX_' extractfield]),...
                    dataIN(ippant,iqntl).(['highX_' extractfield])];
                fullY = [dataIN(ippant,iqntl).(['lowY_' extractfield]), dataIN(ippant,iqntl).(['mainY_' extractfield]),...
                    dataIN(ippant,iqntl).(['highY_' extractfield])];

                %%>> Store
            GFX_fitsObserved_perQntl(ippant,iqntl,1,:)= fullX;
            GFX_fitsObserved_perQntl(ippant,iqntl,2,:)= fullY;
            GFX_stimlists(ippant,iqntl,:) = StimList;
            GFX_params_qntl(ippant, iqntl,:)= [dt(1), dt(2), dt(6)]; % thresh and slope (thresh at 50%).
            GFX_Acc_perQntl(ippant, iqntl,:)= mAcc;
            end % per qntl
        end % ippant.
        


    %% outlier removal.
    if cleanOnly==1
    [ppOutliers_slp,ppOutliers_thrsh]=deal([]); 
    % focusing on the slope and thresh.
    for iqntl=1:5
                        ppOutliers_thrsh(:,iqntl)= isoutlier(GFX_params_qntl(:, iqntl,1));
                        ppOutliers_slp(:,iqntl)= isoutlier(GFX_params_qntl(:,iqntl,2));
            
    end
    % how many in total per pp?
    ppOutCount = sum(ppOutliers_slp,2) +sum(ppOutliers_thrsh,2) ;


    badppant = ppOutCount>0;

useppants= 1:size(GFX_Acc_perQntl,1);
useppants(badppant)=[];

% aggressive rejection:
plist= [1 2 3 4 5 6 7 8 10 11 12 17 19 21 24 25 27 28 31 33 34];
useppants= plist;
    else
        useppants= 1:size(GFX_Acc_perQntl,1);
    end

%     



        %% % continue with plots: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
subplot(2,3,1); cla
            mP= nanmean(GFX_headdata(useppants,:),1);
%             errP = CousineauSEM(GFX_headdata);
            plot(1:size(GFX_headdata,2),mP, ['k-o'])
            %or fancy version:
            for ip= 1:length(mP);

                plot(ip, mP(ip), 'o-','color',qntlCols(ip,:) );
                hold on
            end
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
                ph = patch('XData', xv, 'YData', yv, 'FaceColor', thisColr, 'FaceAlpha', .2 , 'LineStyle','none');
                %             ph.FaceColor = thisColr;
            end
            ylim(yls);
           
            set(gca, 'fontsize', cfg.fontsize, 'ytick', []);
            ylabel('head height')
            xlabel('% step completion')
           
            
            %% now plot the psychometric curves for each portion.

            jitwth = .02;
            % cla
            legh=[];
            % overall Stim list: % if fitting to group data.
            GStimList = squeeze(mean(mean(GFX_stimlists(useppants,:,:),1),2));

            allGParams=[];
            offsets = [-.002:.001:.002];
            for iqntl= 1:5
                figure(1)
                thisColr= qntlCols(qntlBounds(iqntl),:);
                %per quantile, add a small amount of jitter.
                %                 jitt = (rand(1)-.5)/10;
                subplot(2,3,[2,3])
                hold on;

                StimList= squeeze(mean(GFX_stimlists(useppants,iqntl,:),1));
                %revert to full search.
                ProportionCorrectObserved = squeeze(nanmean(GFX_Acc_perQntl(useppants, iqntl,:),1));

                hold on
                errObs = CousineauSEM(squeeze(GFX_Acc_perQntl(useppants,iqntl,:)));
               

                if applyOffset
                    StimList= StimList+offsets(iqntl);
                    GStimList= GStimList+offsets(iqntl);
                end
                if FitatGroup==1 || shareX==1
                legh(iqntl)=plot(GStimList,ProportionCorrectObserved,['.'],'color',thisColr,'markersize',20);
                errorbar(GStimList, ProportionCorrectObserved, errObs, 'k', 'linestyle', 'none');

                else
                legh(iqntl)=plot(StimList,ProportionCorrectObserved,['.'],'color',thisColr,'markersize',20);
                errorbar(StimList, ProportionCorrectObserved, errObs, 'k', 'linestyle', 'none');
    
                end
                
                % add a group fit (mean of individual fits):
                av_IndFit_X = squeeze(mean(GFX_fitsObserved_perQntl(useppants,iqntl,1,:),1));
                av_IndFit = squeeze(mean(GFX_fitsObserved_perQntl(useppants,iqntl,2,:),1));
                err_IndFit = CousineauSEM(squeeze(GFX_fitsObserved_perQntl(useppants,iqntl,2,:)));


                if FitatGroup==1 % compare to high group level fits.

                    %% fit options:
                    options =[];
                    options.sigmoidName = 'norm';   % choose a cumulative Gaussian as the sigmoid
                    options.expType= 'YesNo'; % free params for lapse and guess rate.

                    dataFit = [GStimList, ProportionCorrectObserved, repmat(100,[7,1])];

                    % perform fit:
                    fitresult = psignifit(dataFit,options);
                    resultSmall = rmfield(fitresult,{'Posterior','weight'});
                    clear fitresult % keep mem free as we go
                    % fake plot!
                    set(gcf,'visible', 'off'); % hide figure, but we want the output:
                    figure(4)
                    pr=plotPsych(resultSmall);
                    pGroupX = [pr.low.XData,pr.main.XData,pr.high.XData];
                    pGroupY = [pr.low.YData,pr.main.YData,pr.high.YData];

                    figure(1);  subplot(2,3,[2,3])
                    if applyOffset ==1
                        pGroupX = pGroupX+offsets(iqntl);
                    end
                    plot(pGroupX, pGroupY, '-', 'color',thisColr , 'linew', 2);

                    %result.fit = T W L G;
                    %     slope = (1 - G - L) / (4 * W) **( at threshold)
                    slope = (1 - resultSmall.Fit(4) - resultSmall.Fit(3)) / (4* resultSmall.Fit(2));
                    disp(['Group slope q ' num2str(iqntl) ' ' num2str(slope)]);
                    allGParams(iqntl, 1)= resultSmall.Fit(1); %thresh
                    allGParams(iqntl,2) =resultSmall.Fit(2); % width
                    allGParams(iqntl,3)=slope; % slope

                else

%                     plot(av_IndFit_X, av_IndFit, '-', 'color',thisColr , 'linew', 2);
                    shadedErrorBar(av_IndFit_X, av_IndFit, err_IndFit, { 'linestyle', '-', 'color',thisColr  , 'linew', 1},1);

                end



                set(gca,'color', 'w')
    
                


                if debugPlots==1


                    figure(2);
                    subplot(1,5,iqntl);
                    for ippant= useppants;
                        hold on;
                       pr= plot(squeeze(GFX_fitsObserved_perQntl(ippant,iqntl,1,:)), squeeze(GFX_fitsObserved_perQntl(ippant,iqntl,2,:)));
                       if mod(ippant,2)==0; 
                       text(pr.XData(1), pr.YData(1), ['pp ' num2str(ippant)]);
                       else
                           text(pr.XData(end), pr.YData(end), ['pp ' num2str(ippant)]);

                       end
                    end
                end

            end % iqnt
            %%
%
            %plot error (horizontal error, spread of ppant bins)
%                     errorbar(StimList_Group, repmat(.35, 1, length(StimList_Group)),...
%                         StimList_Group_err, 'k', 'horizontal','linestyle','none', 'LineWidth',1)
% add markers to clarify
%
% width = [-.05*jitwth, .05 *jitwth];
% % hghtMax = [.375 .5 .625 .75 .825 .855 .875];
% for iH = 1:length(StimList_Group)
% %draw vertical line to data.
% vertLims = [.35 hghtMax(iH)];
% plot([StimList_Group(iH) StimList_Group(iH)], vertLims, 'k:')
% 
% end
%
figure(1)
legend(legh, {'qntl1', 'qntl2', 'qntl3', 'qntl4', 'qntl5'}, 'location',' NorthWest', 'autoupdate', 'off');

%%
             %add text with relv stats: (need rmANOVA)
        % convert data to table:
        subjs = [1:length(useppants)]';
        [Fs,Ps]=deal(zeros(2,1));
        
        txtmsgs={};
        for iparam=1:3 % for slope and threshold
            %make columns of the repeated measures (each quantile).
            measures=  squeeze(GFX_params_qntl(useppants, :,iparam));       
             t=splitvars(table(subjs, measures));
             %to make it robust to toggling gait quantiles, extract the
             %autonames:
             %wilkinson notation for a rmmodel.
             rmmodel = ['measures_1-measures_' num2str(length(qntlBounds)-1) '~1'];
        
             WithinDesign = table([1:size(measures,2)]','VariableNames',{'quantiles'});

             mdlfit= fitrm(t, rmmodel, 'WithinDesign', WithinDesign);
             
             %get stats from table:
             rtable = ranova(mdlfit);
             Fs(iparam,1) = rtable.F(1);
             Ps(iparam,1) = rtable.pValue(1);
        textmsgs{iparam} = ['\rm\itF\rm(' num2str(rtable.DF(1)) ',' num2str(rtable.DF(2)) ')=' sprintf('%.3f',Fs(iparam)) ', \itp\rm=' sprintf('%.3f', Ps(iparam))];


        end
           
        text(4,.2,textmsgs, 'fontsize', cfg.fontsize);
        
%             xlim([.5 7.5])
            xlabel('contrast')
           
           
            
        if useLogscale==1
       xlabel('Target contrast (dB)');       
       
        set(gca, 'Xtick',GStimList,'xticklabel', num2str(round(20*log10(GStimList),1)));
        %%
        else
        xlabel('contrast increment');
        set(gca, 'Xtick',GStimList, 'color', [.8 .8 .8]);
            
        end
        ylabel('proportion correct');
        set(gca, 'fontsize', cfg.fontsize);
        
        
        %% add summary bar chart to clarify effect:
        % grouped bar chart:
        barDataThresh = squeeze(GFX_params_qntl(useppants, :,1));
        barDataWidth= squeeze(GFX_params_qntl(useppants, :,2));   
        barDataSlope = squeeze(GFX_params_qntl(useppants, :,3));
        
        barDataPlot = {barDataThresh, barDataWidth, barDataSlope};
        ytitles = {'threshold values',  'width parameter','slope parameter'};
        
        ylimsare=[0,.6; 0,200];
        subspots=[3,0,4];
        
        for iparam= 1:3 
            if iparam==2
                continue % skip the width data
            end
          subplot(2,2 ,subspots(iparam));
          cla
        barD = barDataPlot{iparam};
           mBar = nanmean(barD,1);
       bh= bar(1:size(barD,2),mBar,'FaceColor','flat'); %flat allows us to update colours below
       bh.CData = qntlCols(qntlBounds(1:end-1),:);
       hold on;
       stE = CousineauSEM(barD);
       errorbar(1:size(barD,2),mBar, stE,'linestyle', 'none', 'color', 'k', 'linew',2)
%        ylim(ylimsare(ithrSlp,:))
% adjust ylims to sdrange.
        sdrange = max(mBar) - min(mBar);
        ylim([min(mBar)-1.5*sdrange max(mBar)+1.5*sdrange])

        yts= get(gca, 'ylim');
        p80=yts(1) + .8*(yts(2)-yts(1)); % 80% the height of the subplot.
        %%
        plot([1, 5], [p80,p80], 'k-', ...
            'linew', 2);
        
        if iparam~=1
        text(3, p80, textmsgs{iparam}, 'HorizontalAlignment', 'center',...
            'fontsize', cfg.fontsize,...
            'VerticalAlignment', 'bottom');
        else
            
        text(3, p80, ['\itns'], 'HorizontalAlignment', 'center',...
            'fontsize', cfg.fontsize,...
            'VerticalAlignment', 'bottom');
        end
        %%
if FitatGroup==1
    yyaxis right
    bar(allGParams(:, iparam), 'FaceAlpha',.2)
end
       title(ytitles{iparam});
       ylabel(ytitles{iparam});
       xlabel('quintile');
       set(gca,'fontsize',cfg.fontsize)
%        axis tight
        end % each param to plot.

%%

        cd([cfg.figdir filesep  cfg.type ' onset Qntl PsychFits'])
        
        print([psubj ' ' cfg.type ' onset qntl PsychFits_MS'],'-dpng');
        
end% GFX
end