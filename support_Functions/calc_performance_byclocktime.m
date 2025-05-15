% calc_performance_byclocktime

% in response to a reviewer request, show performance over time (both
% conditions).

%load ppantTables. display target onset distribution and accuracy per
%target location (binned in 100 ms steps).



job=[];
job.concatGFX=0;
job.createShuff=0;
job.testShuff=1;
job.plotResults=1;



setmydirs_detectv3;
% show ppant numbers:
cd(procdatadir)
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%

%load the GFX observed for comparison.
cd('GFX');
load('GFX_Data_inGaits_FourierFits.mat', 'GFX_FourierNull');

%NB after review, now including trial rejection if flagged as outside of
%acceptable DVA (i.e. gaze was elsewhere).
% This is precomputed in AA_TargetLocation_inTrial;

includeShuffle=0;

%%
if job.concatGFX==1
    % preallocate
    [GFX_clocktime_binAcc_stnd,...
        GFX_clocktime_binRT_stnd,...
        GFX_clocktime_binResp_stnd,...
        GFX_clocktime_binTargCounts_stnd,...
        GFX_clocktime_binAcc_wlk,...
        GFX_clocktime_binRT_wlk,...
        GFX_clocktime_binResp_wlk,...
        GFX_clocktime_binTargCounts_wlk] = deal([]);


    for ippant =1:nsubs
        cd(procdatadir)
        load(pfols(ippant).name, 'subjID', 'trial_summaryTable');


        % first store the histcounts of all target onsets, spaced in 100 ms.

        allTargetOnsets = round(trial_summaryTable.targOnset,2);
        allRespOnsets = round(trial_summaryTable.targRT,2);

        bins = 0.01:.01:9.01; % avoid counting zeros (in resp data).

        standTrials = find(trial_summaryTable.isStationary==1);
        walkTrials = find(trial_summaryTable.isStationary==0);

        allTs_wlk = allTargetOnsets(walkTrials);
        allTs_stand = allTargetOnsets(standTrials);

        allRs_wlk = allRespOnsets(walkTrials);
        allRs_stand = allRespOnsets(standTrials);

        bin_Acc_wlk=[];
        bin_RT_wlk=[];
        bin_Resp_wlk=[];


        bin_Acc_stnd=[];
        bin_RT_stnd=[];
        bin_Resp_stnd=[];

        for ibin= 1:length(bins)

            % all target onsets this bin:
            binrow= find(allTargetOnsets==bins(ibin));

            % intersect for performance when standing/walking:
            aa= intersect(standTrials, binrow);
            bb= intersect(walkTrials, binrow);

            bin_RT_stnd(ibin) = nanmean(trial_summaryTable.clickRT(aa));
            bin_Acc_stnd(ibin) = nanmean(trial_summaryTable.targCor(aa));
            bin_counts_stnd(ibin) = length(aa);


            bin_RT_wlk(ibin) = nanmean(trial_summaryTable.clickRT(bb));
            bin_Acc_wlk(ibin) = nanmean(trial_summaryTable.targCor(bb));
            bin_counts_wlk(ibin) = length(bb);


            % all resp onsets this bin:
            binrow= find(allRespOnsets==bins(ibin));
            % intersect for performance when standing/walking:
            aa= intersect(standTrials, binrow);
            bb= intersect(walkTrials, binrow);

            bin_Resp_stnd(ibin)= length(aa);
            bin_Resp_wlk(ibin)= length(bb);

        end

        % average per 20 ms, to approximate main analysis.
        pidx= 1:2:length(bins);
        usedata={bin_Acc_stnd, bin_RT_stnd, bin_Resp_stnd, ...
            bin_Acc_wlk, bin_RT_wlk, bin_Resp_wlk,...
            bin_counts_stnd, bin_counts_wlk};

        for itype= 1:8
            out_bin=nan(1,length(pidx)-1);


            datatobin = usedata{itype};

            for ibin=1:length(pidx)-1
                idx = pidx(ibin):pidx(ibin+1);

                out_bin(ibin) = mean(datatobin(idx), 'omitnan');

            end
            switch itype
                case 1
                    GFX_clocktime_binAcc_stnd(ippant,:)= out_bin;
                case 2
                    GFX_clocktime_binRT_stnd(ippant,:)= out_bin;
                case 3
                    GFX_clocktime_binResp_stnd(ippant,:)= out_bin;

                case 4
                    GFX_clocktime_binAcc_wlk(ippant,:)= out_bin;
                case 5
                    GFX_clocktime_binRT_wlk(ippant,:)= out_bin;
                case 6
                    GFX_clocktime_binResp_wlk(ippant,:)= out_bin;
                case 7
                    GFX_clocktime_binTargCounts_stnd(ippant,:)= out_bin;
                case 8
                    GFX_clocktime_binTargCounts_wlk(ippant,:)= out_bin;

            end
        end % itype

        disp(['fin clock bin for ippant ' num2str(ippant)]);
    end

    cd('GFX')
%%
    save('GFX_clocktime', 'GFX_clocktime_binAcc_stnd',...
        'GFX_clocktime_binRT_stnd','GFX_clocktime_binResp_stnd',...
        'GFX_clocktime_binAcc_wlk','GFX_clocktime_binRT_wlk','GFX_clocktime_binResp_wlk',...
        'GFX_clocktime_binTargCounts_stnd','GFX_clocktime_binTargCounts_wlk');%,'-append');
end
%% %%%%%%%%%%%%%%%%
if job.createShuff==1

    nPerm=1000;

    % same as above, but shuffle in time.
    % preallocate
    [GFX_clocktime_binAcc_stnd_shuff,...
        GFX_clocktime_binRT_stnd_shuff,...
        GFX_clocktime_binResp_stnd_shuff,...
        GFX_clocktime_binTargCounts_stnd_shuff,...
        GFX_clocktime_binAcc_wlk_shuff,...
        GFX_clocktime_binRT_wlk_shuff,...
        GFX_clocktime_binResp_wlk_shuff,...
        GFX_clocktime_binTargCounts_wlk_shuff] = deal(zeros(36,nPerm, 450));


        %%
        binindex = 1:size(GFX_clocktime_binAcc_stnd,2);
        % we actually need to restrict it to within target onsets.
        tmpT = squeeze(mean(GFX_clocktime_binTargCounts_wlk,1));
        shuffRange = [find(tmpT,1,'first'),find(tmpT,1,'last')];
        for iperm=1:nPerm

            % shuffle in time all data points.
            randindx = shuffle([shuffRange(1):shuffRange(2)-1]);
            newvec = [1:shuffRange(1)-1, randindx, shuffRange(2):length(binindex)]; 

                        GFX_clocktime_binAcc_stnd_shuff(:,iperm,:)= GFX_clocktime_binAcc_stnd(:,newvec);
                        GFX_clocktime_binRT_stnd_shuff(:,iperm,:)= GFX_clocktime_binRT_stnd(:,newvec);
                        GFX_clocktime_binResp_stnd_shuff(:,iperm,:)= GFX_clocktime_binResp_stnd(:, newvec);

                        GFX_clocktime_binAcc_wlk_shuff(:,iperm,:)= GFX_clocktime_binAcc_wlk(:, newvec);
                        GFX_clocktime_binRT_wlk_shuff(:,iperm,:)= GFX_clocktime_binRT_wlk(:, newvec);

                        GFX_clocktime_binResp_wlk_shuff(:,iperm,:)= GFX_clocktime_binResp_wlk(:, newvec);
                        GFX_clocktime_binTargCounts_stnd_shuff(:,iperm,:)= GFX_clocktime_binTargCounts_stnd(:,newvec);
                        GFX_clocktime_binTargCounts_wlk_shuff(:,iperm,:)= GFX_clocktime_binTargCounts_wlk(:,newvec);

                
            if mod(iperm,50)==0
                disp(['fin ' num2str(iperm)  ' perms'])
            end


        end % iperm
        %%
        cd(procdatadir);
cd('GFX');
    save('GFX_clocktime', 'GFX_clocktime_binAcc_stnd_shuff',...
        'GFX_clocktime_binRT_stnd_shuff','GFX_clocktime_binResp_stnd_shuff',...
        'GFX_clocktime_binAcc_wlk_shuff','GFX_clocktime_binRT_wlk_shuff','GFX_clocktime_binResp_wlk_shuff',...
        'GFX_clocktime_binTargCounts_stnd_shuff','GFX_clocktime_binTargCounts_wlk_shuff','-append');
end

if job.testShuff==1

% per perm, bin, then fit fourier, retain CV.

GFX_clocktime_Fourierfits=[];

    useData= {GFX_clocktime_binTargCounts_wlk_shuff,GFX_clocktime_binAcc_wlk_shuff, GFX_clocktime_binRT_wlk_shuff, GFX_clocktime_binResp_wlk_shuff,...
        GFX_clocktime_binTargCounts_stnd_shuff,GFX_clocktime_binAcc_stnd_shuff, GFX_clocktime_binRT_stnd_shuff, GFX_clocktime_binResp_stnd_shuff};


     for itype = 1:length(useData)
    %
    testData = useData{itype};

        %% next calculate the participant average over 1 sec intervals.
        hold on;
        % show overlapping 1 second segments:
        GFX_epoched=[];
        disp('epoching data');
        for ippant = 1:size(testData,1);
            newD=[];
            for ibin= 1:9

                binindex= [1:50] + 50*(ibin-1);

                newD(ibin,:,:)= testData(ippant,:,binindex);

            end
            GFX_epoched(ippant,:,:)= nanmean(newD,1);
        end
        
        %% now perform test at all freqs:

        FitType = 'fourier1';
        % Creating and showing a table array to specify bounds.
        CoeffNames = coeffnames(fittype(FitType));

        %set bounds for w
        CoeffBounds = array2table([-Inf(1,length(CoeffNames));...
            Inf(1,length(CoeffNames))],'RowNames',...
            ["lower bound", "upper bound"],'VariableNames',CoeffNames);
        %
        % Specifying bounds according to the position shown by the table.
        % e.g. to force fit with w ~ 1.545, we ChangeBound of w parameter
        % CoeffBounds.w(1) = 1.54;
        % CoeffBounds.w(2) = 1.55;
        %
        Hzspace = [0.01:.2:12];
        % perW=per
       
        fits_Rsquared_obsrvd = nan(1, length(Hzspace));
        fits_Rsquared_shuff=  nan(1, length(Hzspace));
        groupShuff = squeeze(nanmean(GFX_epoched,1));
        %
        % step through w, forcing fit at particular periods, by updating the
        % bounds in fit options.
        
        xvec= [1:50];

        for ifreq= 1:length(Hzspace)
            % include period and Rsquared
            %treat max xvec as our full 'period'
            %             Hzapp = xvec(end)/ (2*pi/(f.w));
            testw = 2*pi*Hzspace(ifreq)/xvec(end);

            CoeffBounds.w(1) = testw;
            CoeffBounds.w(2) = testw;

            %update fit opts settings

            %Update Fit Options setting.
            FitOpts = fitoptions('Method','NonlinearLeastSquares','Lower',table2array(CoeffBounds(1,:)),...
                'Upper',table2array(CoeffBounds(2,:)));

            %first test this period on observed data
            %set last coefficient value in the fourier model (w) to this value:
        for iperm=1:size(groupShuff,2)
          

            try
            [tmpf,gof] = fit(xvec', squeeze(groupShuff(iperm,:))', 'fourier1', FitOpts);
            fits_Rsquared_shuff(iperm,ifreq) = gof.rsquare;
            catch
                disp([' skipping perm ' num2str(iperm)  ', epochs, likely a NaN'])
            end
        end % per perm
        disp(['fin perms for freq ' num2str(ifreq)]);
        end % per freq

            % now 
  maxperPerm = max(fits_Rsquared_shuff,[],2);
%         fits_Rsquared_shuffCV(:,ifreq) = quantile(fits_Rsquared_shuff(:,ifreq), [.05, .5, .95]);
       fits_Rsquared_shuffCV_new = quantile(maxperPerm, [.05, .5, .95]);
% store the 95% of all freqs:
GFX_clocktime_Fourierfits(itype).fits_Rsquared_shuff = fits_Rsquared_shuffCV_new;
disp(['saving']);
save('GFX_clocktime', 'GFX_clocktime_Fourierfits', '-append');
     end % per type
      
end % job.



%% %%%%%%%%%%%%%%%%%%%%
if job.plotResults==1
    %% plot grand mean result:
    clf

    ylabsAre={'Target onsets', 'Accuracy', 'RT', 'Response counts'};
    DVsare= {'counts','Acc', 'rts','counts'}; % for accessing the prev data.
    Onsetsare= {'Target','Target', 'Target','Response'}; % for accessing the prev data.

    colsAre= {'k','b', 'r', 'm'};
    useData= {GFX_clocktime_binTargCounts_wlk,GFX_clocktime_binAcc_wlk, GFX_clocktime_binRT_wlk, GFX_clocktime_binResp_wlk,...
    GFX_clocktime_binTargCounts_stnd,GFX_clocktime_binAcc_stnd, GFX_clocktime_binRT_stnd, GFX_clocktime_binResp_stnd};

    % subplot(311);
    % bar(bins, nanmean(GFX_clocktime_binTargOns,1));
    % ylabel('Target onset (count)');

for ifg=1:2
    figure(ifg); clf;
    set(gcf,'color','w','units','normalized','position',[0 0 .95 .9])
end

    fntsize=14
    iDVcount=1;
    for iDV=1:4%length(useData)
        if iDV<5 
            figure(1);
            iDVcount=iDV;
        else 
            figure(2);

            iDVcount= iDV-4;
        end

        %% first plot the whole trial DV ( 9 sec).
        subplot(4,4,[1,2]+ 4*(iDVcount-1));cla;

        plotData = useData{iDV};
        bh=bar(pidx(1:end-1), nanmean(plotData,1));
        bh.FaceColor= [.9 .9 .9];
        bh.EdgeColor= colsAre{iDVcount};
        ylabel(ylabsAre{iDVcount})

        if iDVcount==1
            if iDV==1
            title('When Walking');
            else
            title('When Standing');
            end
        end

        if iDV==4
            xlabel('Trial time (sec)');
        end
        set(gca,'xtick', [1 100,200,300,400,500,600,700,800,900],...
            'xticklabel', {'0', '1','2','3','4','5','6','7','8','9'}, ...
            'fontsize',fntsize);

        %% next calculate the participant average over 1 sec intervals.
        hold on;
        % show overlapping 1 second segments:
        GFX_epoched=[];
        for ippant = 1:size(plotData,1);
            newD=[];
            for ibin= 1:9

                binindex= [1:50] + 50*(ibin-1);

                newD(ibin,:)= plotData(ippant,binindex);

            end
            GFX_epoched(ippant,:)= nanmean(newD,1);
        end
        %
        subplot(4,4,3+ 4*(iDVcount-1));cla

        gM=  nanmean(GFX_epoched,1);
        bar(1:50, gM, 'FaceColor',[.9 .9 .9],  'EdgeColor', colsAre{iDVcount}, 'linew',1);
        stE= CousineauSEM(GFX_epoched);

        hold on
        errorbar(1:50, gM, stE, 'color',[.5 .5 .5],'linestyle', 'none', 'linew',2);
        sdrange = max(gM) - min(gM);
        ylim([min(gM)-.5*sdrange max(gM)+1*sdrange])

        ylsat = get(gca, 'ylim');
        if ylsat(1)<0
            ylim([0 ylsat(2)])
        end
        xvec= [1:50];

        % apply best fit (overlay)
        [f,gof]= fit(xvec',gM',  'fourier1');

        % plot:
        hold on;
        %             yyaxis right
        h=plot(f, xvec, gM);%,
        h(2).LineWidth = 2;
        h(2).Color = colsAre{iDVcount};
        %
        %treat max xvec as our full 'period'
        fitperiod = f.w;
        %convert to period per samples.
        % include period and Rsquared
        %treat max xvec as our full 'period'
        Hzapp = xvec(end)/ (2*pi/(f.w));

        % critical value for sig:
         shuffCV=GFX_clocktime_Fourierfits(iDV).fits_Rsquared_shuff(3);

        if gof.rsquare>shuffCV
        legdetails = [sprintf('%.2f', Hzapp) ' Hz, R^2 = ' sprintf('%.2f', gof.rsquare) ];
        else
        legdetails = ['\itns'];
%         clear h
set(h(2),'visible','off')
% set(h(1),'visible','off')
        end
        legend(h(2), legdetails, 'fontsize', 12, 'autoupdate', 'off', 'Location', 'North')

        ylabel(ylabsAre{iDVcount})
        set(gca,'fontsize', fntsize, 'xtick', [1,   25,  xvec(end)], 'XTickLabels', {'0', '.5', '1'})
        xlabel('')
        if iDVcount==4
            xlabel([ 'Event onset (seconds)']);%
        end


        if iDVcount==1
            title('Average (consecutive epochs)');
        end
        %% now perform test at all freqs:

        FitType = 'fourier1';
        % Creating and showing a table array to specify bounds.
        CoeffNames = coeffnames(fittype(FitType));

        %set bounds for w
        CoeffBounds = array2table([-Inf(1,length(CoeffNames));...
            Inf(1,length(CoeffNames))],'RowNames',...
            ["lower bound", "upper bound"],'VariableNames',CoeffNames);
        %
        % Specifying bounds according to the position shown by the table.
        % e.g. to force fit with w ~ 1.545, we ChangeBound of w parameter
        % CoeffBounds.w(1) = 1.54;
        % CoeffBounds.w(2) = 1.55;
        %
        Hzspace = [0.01:.2:12];
        % perW=per
        fits_Rsquared_obsrvd = nan(1, length(Hzspace));

        %best fit (use model params below):
        f = fit(xvec', gM', 'fourier1'); %unbounded
        %
        % step through w, forcing fit at particular periods, by updating the
        % bounds in fit options.
        for ifreq= 1:length(Hzspace)
            % include period and Rsquared
            %treat max xvec as our full 'period'
            %             Hzapp = xvec(end)/ (2*pi/(f.w));
            testw = 2*pi*Hzspace(ifreq)/xvec(end);

            CoeffBounds.w(1) = testw;
            CoeffBounds.w(2) = testw;

            %update fit opts settings

            %Update Fit Options setting.
            FitOpts = fitoptions('Method','NonlinearLeastSquares','Lower',table2array(CoeffBounds(1,:)),...
                'Upper',table2array(CoeffBounds(2,:)));

            %first test this period on observed data
            %set last coefficient value in the fourier model (w) to this value:

            [tmpf,gof] = fit(xvec', gM', 'fourier1', FitOpts);
            % how good/bad was the fit?
            fits_Rsquared_obsrvd(1,ifreq) = gof.rsquare;
        end
        %
        subplot(4,4,4+ 4*(iDVcount-1)); cla
        hzf=plot(Hzspace, fits_Rsquared_obsrvd, colsAre{iDVcount}, 'linew',2);
hold on;
        % add the CV:
        shuffCV=GFX_clocktime_Fourierfits(iDV).fits_Rsquared_shuff(3);
        hzN=plot(xlim,[shuffCV, shuffCV], 'k:', 'linew',2);


        ylim([0 .75])
        ylabel('R^2');
        if iDVcount==4
            xlabel('Frequency (Hz or cps)');
        end
        hold on; % now plot the observed from stride-cycle data.
        % DVsare= {'Acc', 'rt','counts'}; % for accessing the prev data.
        % Onsetsare= {'Target', 'Target','Response'}; % for accessing the prev data.

        compareObs = GFX_FourierNull.([Onsetsare{iDVcount} 'Ons_doubgc_binned_' DVsare{iDVcount} '_fitsRsq_Obs']);

        Hzspace = [0.01:.2:10];
        cps= plot(Hzspace, compareObs, 'color',[.7 .7 .7],'linew',2);
        plot(Hzspace, compareObs, 'color',[.7 .7 .7], 'linew',2);
        set(gca,'fontsize',fntsize);
ylim([0 .75])
        lg=legend([hzf, cps,hzN], {'Fit strength (Hz)', 'Fit strength (cps)', '95% CI'}, 'fontsize', fntsize-4);
%         text(6, .65, 'Fit strength (Hz)', 'fontsize', fntsize-2,'color',colsAre{iDVcount})
            
        if iDVcount==1
            title('Fit strength');
        end
    end
end
%%




