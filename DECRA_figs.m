%DECRA figs

job= [];
job.plotFig1 = 1;  % this is the 3D trial overview.
job.plotFig1_v2 = 0;  % this is the 3D trial overview stripped back (KISS)
job.plotFig2 = 0; % this is the GFX summary, includes sensitivity and aud pilot
job.plotFig2_ver2 = 0; % this is the GFX summary, includes sensitivity and aud pilot, smaller version.
job.plotFig3= 0; % plots the subitizing pilot results (both ranges), single panel.
job.plotFig4_EEG_mobi=0;
%plot_example3Dtrial
job.plotACNSpanel1= 0; % the detectv3 summary, new spacing
% example walk?

clf
setmydirs_detectv3;
cd(procdatadir);
cd('GFX')
%
load('GFX_Data_inGaits.mat','GFX_headY', 'GFX_TargPosData', 'pidx2');
targCounts = [];
targAcc = [];
targRT=[];
for ippant = 1:size(GFX_TargPosData,1)
    targCounts(ippant,:)= GFX_TargPosData(ippant,3).doubgc_binned_counts;
    targAcc(ippant,:)= GFX_TargPosData(ippant,3).doubgc_binned_Acc;
    targRT(ippant,:) = GFX_TargPosData(ippant,3).doubgc_binned_rts;

end

%
headY=[];
for ippant = 1:36
    headY(ippant,:) = GFX_headY(ippant).doubgc;
end
fntsize=16;
%%
if job.plotFig1==1
    cd(procdatadir);
    subjID= 'BL';

    % itrial = 23, 60,84
    itrial=84;
    cd(procdatadir)
    lfile  = dir([pwd filesep subjID '*']);
    load(lfile.name)
    % fntsize= cfg.fntsize;
    
    HeadPos(itrial).Z= HeadPos(itrial).Z- mean(HeadPos(itrial).Z);
    timevec= trialInfo(itrial).times;
    % timevec= HeadPos(itrial).X +5;
    figure(1); clf
%     set(gcf,'units','normalized','position', [0 0 .5 .4]); 
    set(gcf,'units','normalized','position', [0 0 .56 .4]); 
    
    figure(2); clf
    set(gcf,'units','normalized','position', [0 0 .4 .4]);
    %
    figure(1);
    hold on;
    p3D=plot3(nan, nan,nan, ['m-'], 'LineWidth',8);

    plot3(timevec, HeadPos(itrial).Z,HeadPos(itrial).Y, ['m-'], 'LineWidth',4);

    %trick to add shading.
    x= timevec;
    % y= HeadPos(itrial).Z;
    % z= HeadPos(itrial).Y;
    % r = sqrt(x.^2 + y.^2 + z.^2);
    % g = patch('Vertices', [x(:), y(:),z(:); nan nan nan], 'Faces', (1:length(x)+1).', 'FaceVertexCData', [r(:); nan], 'EdgeColor', 'm', 'Marker','.','linew',2);
    %
    % p.Alpha=.2
    % axis tight;
    backwall = .15;
    floorwall = 1.6; % 1.6
    xlabel('Trial time (sec)');
    % ylabel('Head sway (m)', 'rotation', -35);
    set(gca,'ytick',[])
    ylim([-backwall backwall])
    zlabel('Head height (m)');
    zlim([floorwall floorwall+.15])
    %
    set(gca,'fontsize', fntsize);
    hold on;
    %plot height:
    plot3(timevec, repmat(backwall, [1 length(trialInfo(itrial).times)]), HeadPos(itrial).Y, 'color', [.7 .7 .7], 'linew', 2);
    % ylim([1.69 1.78])
    hold on;
    % nan for legend (plot floor):
    p2D=plot3(timevec, HeadPos(itrial).Z, repmat(floorwall, [1 length(trialInfo(itrial).times)]), 'color', [.7 .7 .7], 'linew', 2);
   %


    grid on
    view([-15 38]); %
    view([-5.2 54])
    view([-10.8, 23.45])
    view([-38.2, 29])
    set(gca,'fontsize',fntsize, 'xtick', [0:1:10], 'xticklabels', ...
        {'0',...'',...
        '1',...'',...
        '2',...'',...
        '3',...''...
        '4',...''...
        '5',...''...
        '6',...''...
        '7',...''...
        '8',...''...
        '9',...''...
        '10'})
    %

    % add markers.

    hold on;

    % add target results?
    %
    % tOnsets = HeadPos(124).tr
    relT= find(trial_summaryTable.trial == itrial);
    tOnsets = trial_summaryTable.targOnset(relT);
    tResult = trial_summaryTable.targCor(relT);
    tResult(4)=0;
    hold on;
    trialTimes = trialInfo(itrial).times;
    HeadData= HeadPos(itrial).Y;
    txtHeight = 1.7;
    for itarg = 1:length(tOnsets)

        pAt = dsearchn(trialTimes, tOnsets(itarg)');

        if tResult(itarg)==1% Hit.
            fillcol=  [.7 .7 .7];
            fillcol3=  'b';
            fillcolx='k';
        else
            fillcol='w';
            fillcol3='w';
            fillcolx='w';
        end

        %add to back wall:
        tl=plot3(timevec(pAt), backwall, HeadData(pAt), 'o', 'color', [.7 .7 .7], 'MarkerFaceColor', fillcol,'linew',2, 'markersize', 10);
        % %and floor:
        plot3(timevec(pAt), HeadPos(itrial).Z(pAt), floorwall, 'o','color', [.7 .7 .7],'MarkerFaceColor',fillcol, 'linew',2, 'markersize', 10);

        %
        %and 3D trace?
        tl=    plot3(timevec(pAt), HeadPos(itrial).Z(pAt), HeadPos(itrial).Y(pAt), 'o', 'color', 'b', 'MarkerFaceColor',fillcol3,'linew',2,'markersize',10);

        if tResult(itarg)==1 % store legend for hit and miss
            hl = tl;
        else
            ml=tl;
        end

        %add to marker on time axis:
        plot3([timevec(pAt)], [ -.15], [floorwall+.005], 'ko', 'linew',2,'MarkerFaceColor', fillcol3,'markersize',10)
        % verical orientation:
        plot3([timevec(pAt),timevec(pAt)], [-.15 -.15], [floorwall,floorwall+.025], 'k-', 'linew',2)

        %
        if itarg==1

            %link to backwall:
            plot3([timevec(pAt) timevec(pAt)], [backwall, HeadPos(itrial).Z(pAt)], [HeadPos(itrial).Y(pAt),HeadPos(itrial).Y(pAt)],['-.'], 'color', ['k'])
            %link to floor:
            plot3([timevec(pAt) timevec(pAt)], [ HeadPos(itrial).Z(pAt), HeadPos(itrial).Z(pAt)], [floorwall,HeadPos(itrial).Y(pAt)],['-.'], 'color', 'k')
        end
    end

    [lgd,hobh]= legend([p3D, hl,ml],{'3D Head position',  'Target Hit', 'Target Miss'},...
        'autoupdate','off', 'Location','North', 'fontsize', 14, 'Position',[0.3378 0.7869 0.2240 0.1991]);

    %shrink size of line obks
    linesh = findobj(hobh, 'type', 'Line');
    % restrict to 50 % size, centred on original centre.
    for ih=[1,3]
        cntr=mean(linesh(ih).XData);
        adj = diff(linesh(ih).XData)/4;
        linesh(ih).XData= [cntr-adj, cntr+adj]
    end
    %  tt=text(0, -.15, 1.65, 'Target onsets:', 'fontsize', 20)
    tt=text(-.25, -.15, 1.63, sprintf('Target \nOnsets:'), 'fontsize', 14, 'Interpreter', 'tex')

    tt=text(-1, backwall, 1.80, 'Walking direction\rightarrow ', 'fontsize', 18,'BackgroundColor','w')
    % tt=text(1, backwall, 1.76, 'Walking direction', 'fontsize', 20,'BackgroundColor','w')
    %draw an arrow?
    % plot3([1 3], [backwall backwall], [1.72 1.72 ], 'k-', 'linew',4)
    % plot3([3 3], [backwall backwall], [1.72 1.72], 'k>', 'linew',4)
    %  text(0.1, 1.71, {['Target'];['onset:']}, 'HorizontalAlignment', 'left', 'fontsize', fntsize-2);
    xlim([0 9])
    box on

    % shadow at floor?
    plot3([timevec(1) timevec(end)], [backwall backwall], [floorwall floorwall], '-', 'color',[.2 .2 .2 .1], 'linew',6)
    plot3([timevec(end) timevec(end)], [backwall backwall], [floorwall floorwall+.15], '-', 'color',[.2 .2 .2 .1], 'linew',6)
    plot3([timevec(end) timevec(end)], [backwall -backwall], [floorwall floorwall], '-', 'color',[.2 .2 .2 .1], 'linew',6)

    %    set(gca,'PlotBoxAspectRatio', [1,.3 .4])
    %% print transparent version:
    set(gca,'color','w'); set(gcf,'color','w'); % 'none' to output a transparent figure.
    print(gcf, '-dpdf', 'fig1c')
%%
    print(gcf, '-dpdf', 'fig1c-pdf', '-painters')
%%
    export_fig 'NComms1A.png' -transparent

    export_fig 'tmppdf.pdf'; 

%% another try:
exportgraphics(gcf, 'fig1c.pdf', 'ContentType','vector')

%     export_fig 'NComms1c.pdf'  % pdf is output with strange errors. 
% export_fig NComms1c.eps
    %
    %%  now plot the head height change
    figure(2);

    subplot(2,2,[1,3]); cla

    % yyaxis right;
    % tmp = nanmean(targCounts,1);
    bh=bar(1:40, nanmean(targCounts,1));
    bh.FaceColor= [.9 .9 .9];
    bh.EdgeColor=[.5 .5 .5]; hold on
    bh.FaceAlpha = .2;
    ylim([0 20])
    set(gca,'YColor', 'k','fontsize',12);
    ylabel('Target counts');%, 'HorizontalAlignment','right')
    % yyaxis left
    plot([1:40], imresize(mean(headY,1),[1,40])*100+20, ['ko-'],'MarkerSize',10,'MarkerIndices',1:1:200);
    ylim([-4 25])
    hold on
    ht= plot(nan, nan, ['ko'],'MarkerSize',10,'MarkerIndices',1:2:200)
    xlabel('');% stride-cycle');
    legend(ht, 'Head height','fontsize', fntsize, 'units','normalized','Position', [0.196 .89 .21 .07])
    % ylabel(['Head height (m)'],'fontsize',13,'HorizontalAlignment','left');
    set(gca,'ytick',[0:2:18],'XTickLabel',[],'YColor','k', 'fontsize', fntsize);
    box off
    %
    %
    %
    hold on
    % subplot(2,2,3);
    % plot([1:200]./2, repmat(1,[1,200]), 'w')
    xlabel('stride-cycle (%)');
    % set(gca,'ytick',[], 'YColor', 'w')

    % set(gca,'fontsize',fntsize);
end

if job.plotFig1_v2==1
    
    cd(procdatadir);
    subjID= 'BL';

    % itrial = 23, 60,84
    itrial=84;
    cd(procdatadir)
    lfile  = dir([pwd filesep subjID '*']);
    load(lfile.name)
    % fntsize= cfg.fntsize;
    
    HeadPos(itrial).Z= HeadPos(itrial).Z- mean(HeadPos(itrial).Z);
    timevec= trialInfo(itrial).times;
    % timevec= HeadPos(itrial).X +5;
    figure(1); clf
    set(gcf,'units','normalized','position', [0 0 .5 .4]);
    figure(2); clf
    set(gcf,'units','normalized','position', [0 0 .4 .4]);
    %
    figure(1);
    hold on;
    % plot 3D head?
%     p3D=plot3(nan, nan,nan, ['m-'], 'LineWidth',8);
%     plot3(timevec, HeadPos(itrial).Z,HeadPos(itrial).Y, ['m-'], 'LineWidth',4);
% room specs:   
    backwall = .15;
    floorwall = 1.6; % 1.6
    xlabel('Trial time (sec)');
    % ylabel('Head sway (m)', 'rotation', -35);
    set(gca,'ytick',[])
    ylim([-backwall backwall])
    zlabel('Head height (m)');
    zlim([floorwall floorwall+.15])
    %
    set(gca,'fontsize', fntsize);
    hold on;
    %plot height:
    plot3(timevec, repmat(backwall, [1 length(trialInfo(itrial).times)]), HeadPos(itrial).Y, 'color', [.7 .7 .7], 'linew', 2);
    hold on;
    
    %plot floor sway?
%     p2D=plot3(timevec, HeadPos(itrial).Z, repmat(floorwall, [1 length(trialInfo(itrial).times)]), 'color', [.7 .7 .7], 'linew', 2);
   %


    grid on
    view([-15 38]); %
    view([-5.2 54])
    view([-10.8, 23.45])
    view([-38.2, 29])
    set(gca,'fontsize',fntsize, 'xtick', [0:1:10], 'xticklabels', ...
        {'0','1','2','3','4','5','6','7','8','9','10'})
    shg
    %

    % add markers.

    hold on;
%     plot3([timevec(83) timevec(83)], [.05, HeadPos(itrial).Z(83)], [HeadPos(itrial).Y(83),HeadPos(itrial).Y(83)] )

    % add target results?
    %
    % tOnsets = HeadPos(124).tr
    relT= find(trial_summaryTable.trial == itrial);
    tOnsets = trial_summaryTable.targOnset(relT);
    tResult = trial_summaryTable.targCor(relT);
    tResult(4)=0;
    hold on;
    trialTimes = trialInfo(itrial).times;
    HeadData= HeadPos(itrial).Y;
    txtHeight = 1.7;
    for itarg = 1:length(tOnsets)

        pAt = dsearchn(trialTimes, tOnsets(itarg)');

        if tResult(itarg)==1% Hit.
            fillcol=  [.7 .7 .7];
            fillcol3=  'b';
            fillcolx='k';
        else
            fillcol='w';
            fillcol3='w';
            fillcolx='w';
        end

        %add to back wall:
        tl=plot3(timevec(pAt), backwall, HeadData(pAt), 'o', 'color', [.7 .7 .7], 'MarkerFaceColor', fillcol,'linew',2, 'markersize', 10);
        % %and floor:
%         plot3(timevec(pAt), HeadPos(itrial).Z(pAt), floorwall, 'o','color', [.7 .7 .7],'MarkerFaceColor',fillcol, 'linew',2, 'markersize', 10);

        %
        %and to 3D trace?
%         tl=    plot3(timevec(pAt), HeadPos(itrial).Z(pAt), HeadPos(itrial).Y(pAt), 'o', 'color', 'b', 'MarkerFaceColor',fillcol3,'linew',2,'markersize',10);
% 
        
        %add to marker on time axis:
        tl=plot3([timevec(pAt)], [ -.15], [floorwall+.005], 'ko', 'linew',2,'MarkerFaceColor', fillcol3,'markersize',10);
        % verical orientation:
        plot3([timevec(pAt),timevec(pAt)], [-.15 -.15], [floorwall,floorwall+.025], 'k-', 'linew',2)
        if tResult(itarg)==1 % store legend for hit and miss
            hl = tl;
        else
            ml=tl;
        end
        
    end

    %legend (3D elements):
%     [lgd,hobh]= legend([p3D, hl,ml], {'3D Head position',  'Target Hit', 'Target Miss'}, 'autoupdate','off', 'Location','North', 'fontsize', 16, 'Position',[0.3378 0.7869 0.2240 0.1991]);
    %legend (2D elements):
    [lgd,hobh]= legend([ hl,ml], {'Target Hit', 'Target Miss'}, 'autoupdate','off', 'Location','North', 'fontsize', 16, 'Position',[0.3378 0.7869 0.2240 0.1991]);

%     %shrink size of line obks
%     linesh = findobj(hobh, 'type', 'Line');
%     % restrict to 50 % size, centred on original centre.
%     for ih=[1,3]
%         cntr=mean(linesh(ih).XData);
%         adj = diff(linesh(ih).XData)/4;
%         linesh(ih).XData= [cntr-adj, cntr+adj]
%     end
    %  tt=text(0, -.15, 1.65, 'Target onsets:', 'fontsize', 20)

    tt=text(-.25, -.15, 1.63, sprintf('Target \nOnsets:'), 'fontsize', fntsize, 'Interpreter', 'tex')
    % tt=text(1, backwall, 1.76, 'Walking direction', 'fontsize', 20,'BackgroundColor','w')
    %draw an arrow on backwall:
    % plot3([1 3], [backwall backwall], [1.72 1.72 ], 'k-', 'linew',4)
    % plot3([3 3], [backwall backwall], [1.72 1.72], 'k>', 'linew',4)
    %  text(0.1, 1.71, {['Target'];['onset:']}, 'HorizontalAlignment', 'left', 'fontsize', fntsize-2);
    
    %draw arrow on floor:
    plot3([1 8], [0 0], [floorwall floorwall ], 'k-', 'linew',4)
    plot3([8 8], [0 0], [floorwall floorwall], 'k>', 'linew',4)
    tt=text(-1, backwall, 1.80, 'Walking direction\rightarrow ', 'fontsize', 20,'BackgroundColor','w')

%     tt=text(-1, 0, floorwall, 'Walking direction\rightarrow ', 'fontsize', 20,'BackgroundColor','w');
    

    xlim([0 9])
    box on

    % shadow at floor?
    plot3([timevec(1) timevec(end)], [backwall backwall], [floorwall floorwall], '-', 'color',[.2 .2 .2 .1], 'linew',6)
    plot3([timevec(end) timevec(end)], [backwall backwall], [floorwall floorwall+.15], '-', 'color',[.2 .2 .2 .1], 'linew',6)
    plot3([timevec(end) timevec(end)], [backwall -backwall], [floorwall floorwall], '-', 'color',[.2 .2 .2 .1], 'linew',6)

    %    set(gca,'PlotBoxAspectRatio', [1,.3 .4])
    %% print transparent version:
    set(gca,'color','none'); set(gcf,'color','none');
    export_fig 'DECRA1av_2.png' -transparent
    end
    %

if job.plotFig2==1
    %% this job also needs the auditory sensitivity data loaded!

    setmydirs_detectv3;
    cd(procdatadir);
    cd(procdatadir); cd('GFX')
    load('GFX_Data_inGaits_FourierFits.mat');
    %%
    figure(1); clf; set(gcf,'color','w')
    % plot variation in accuracy and rt
    % for ippant = 1:length(PFX_FourierNull);
    %     subplot(6,6,ippant);
    %     plot(PFX_FourierNull(ippant).TargetOns_doubgc_binned_rts_fitsRsq_Obs, 'r'); hold on
    %         plot(PFX_FourierNull(ippant).TargetOns_doubgc_binned_Acc_fitsRsq_Obs, 'b'); hold on
    %         title(num2str(ippant))
    % end
    %

    %% first plot the group result
    figure(2); clf
    set(gcf,'units','normalized','position', [.1 .1 .6 .4]);
    fntsize=12;
    useD = {targAcc*100, targRT};
    bars ={'Accuracy', 'Reaction times'};
    ys = {'(%)', '(sec)'};
    cols = {'b', 'r'};
    subsp= [1,5];

    for idata=1:2
        subplot(2,4,subsp(idata))%[4,8])
        % accuracy?
        mdiff = round(mean(diff(pidx2)./2));
        xvec = pidx2(1:end-1) + mdiff;

        tmpD= useD{idata};
        % gM= nanmean(targAcc,1);
        [bh,Fh]=bar_wFourier(xvec, tmpD);
        bh.FaceColor= [.9 .9 .9];
        bh.EdgeColor=cols{idata}; hold on
        bh.FaceAlpha = .2;
        gM = nanmean(tmpD);
        sdrange = max(gM) - min(gM);
        ylim([min(gM)-.5*sdrange max(gM)+.5*sdrange])
        ylsat = get(gca, 'ylim');
        legend([bh Fh(2)], {'data', 'fit'}, 'fontsize',fntsize-4);
        Fh(2).LineWidth=4;
        Fh(2).Color= cols{idata}
        ylabel(ys{idata}, 'units','normalized', 'Position',[-.1 1 1], 'Rotation', 0);%VerticalAlignment','top')
        ya= get(gca,'ylim');
        % place text at 95% mark.
        yAt = (ya(2) - ya(1))*1.1 + ya(1);
        text(5, yAt, bars{idata}, 'fontsize',fntsize)
        % text(bars{idata}, 'HorizontalAlignment','right')
        if idata==2
            xlabel('stride-cycle (%)');
        else
            xlabel(' ')
        end
        set(gca,'fontsize',fntsize)
        %
        hghts= [.875 , .4];
        % legend([Fh(2)], {'group fit'}, 'fontsize',fntsize-4, 'NumColumns',2, 'units','normalized', 'Position', [0.175, hghts(idata),.1 .05]);%,'Position', [0.2 0.2 .1 .1]);
        legend([Fh(2)], {'group fit'}, 'fontsize',fntsize-4, 'NumColumns',2, 'units','normalized', 'Location', 'NorthEast')%'Position', [0.175, hghts(idata),.1 .05]);%,'Position', [0.2 0.2 .1 .1]);

        shg

    end
    set(gcf,'color','w')

    %who some examples of ppant fits.
    % ppant 6 vs 29 for 2 vs 4 cps in Accuracy;
    % clf
    accpp = [20,29];

    % offset =

    % cla
    hghts= [.45, .92];
    for ipp=1:2
        subplot(2,4,2); hold on
        pdata = targAcc(accpp(ipp),:);
        % offset
        pdata= (pdata-mean(pdata)) + ipp/2;
        plot(xvec, pdata, 'ko-');
        % add fit:

        [f,gof]= fit(xvec',pdata',  'fourier1');

        % plot:
        hold on;
        %             yyaxis right
        Fh=plot(f, xvec, pdata);%,
        Fh(2).LineWidth = 2;
        Fh(2).Color = 'b';
        Fh(1).Visible= 'off'
        lgh=Fh(2);
        text(1, ipp/2+.25, ['ppt. ' num2str(ipp)], 'fontsize', fntsize-2,'color','b');
        % save p data for next plot

        if ipp==1
            YD = Fh(2).YData;

        end
        yla= get(gca,'ylim'); % extend by 20%
        dY= yla(2) - yla(1);
        ylim([yla(1) yla(1)+dY*1.75])
        % legend(lgh,'participant fit','fontsize',fntsize-4, 'NumColumns',2, 'units','normalized', 'Position', [0.51, hghts(2),.1 .05]);%,'Position', [0.2 0.2 .1 .1]);
        legend(lgh,'participant fits','fontsize',fntsize-4, 'NumColumns',2, 'units','normalized', 'Location', 'North');% [0.51, hghts(2),.1 .05]);%,'Position', [0.2 0.2 .1 .1]);
    end
    set(gca,'YTickLabel',[])
    ylabel('Accuracy (%)')
    set(gca,'fontsize',fntsize)
    shg

    % %% now same for Rt, focus on amp differences.
    % figure(1); clf
    % for ippant = 1:length(PFX_FourierNull);
    %     subplot(6,6,ippant);
    %    pdata = targRT(ippant,:);
    % [f,gof]= fit(xvec',pdata',  'fourier1');
    %
    % % plot:
    % hold on;
    % %             yyaxis right
    % Fh=plot(f, xvec, pdata);%,
    % Fh(2).LineWidth = 2;
    % Fh(2).Color = 'r';
    % Fh(1).Visible= 'off';
    % legend off
    % ylim([.3 .4])
    % title(num2str(ippant))
    % end
    %
    % ppants 4,13
    RTpp = [4,13];

    % offset =

    % cla
    for ipp=1:2
        subplot(2,4,6); hold on
        pdata = targRT(RTpp(ipp),:);
        % offset
        pdata= (pdata-mean(pdata))*2.5 + ipp/5;
        plot(xvec, pdata, 'ko-');
        % add fit:

        [f,gof]= fit(xvec',pdata',  'fourier1');

        % plot:
        hold on;
        %             yyaxis right
        Fh=plot(f, xvec, pdata);%,
        Fh(2).LineWidth = 2;
        Fh(2).Color = 'r';
        Fh(1).Visible= 'off';
        lgh=Fh(2);
        text(1, ipp/5+.1, ['ppt. ' num2str(ipp+2)], 'fontsize', fntsize-2, 'color','r');

        yla= get(gca,'ylim'); % extend by 20%
        dY= yla(2) - yla(1);
        ylim([yla(1) yla(1)+dY*2])

    end
    % legend(lgh,'participant fit','fontsize',fntsize-4, 'NumColumns',2, 'units','normalized', 'Position', [0.51, hghts(1),.1 .05]);%,'Position', [0.2 0.2 .1 .1]);

    legend(lgh,'participant fits','fontsize',fntsize-4, 'NumColumns',2, 'units','normalized', 'Location','NorthEast');% [0.51, hghts(1),.1 .05]);%,'Position', [0.2 0.2 .1 .1]);

    set(gca,'Ytick',[], 'YTickLabel', {'participant 2', 'participant 1'});
    ylabel('RTs (sec)')
    shg
    shg
    xlabel('stride-cycle (%)');
    set(gca,'fontsize', fntsize)

    %% plot the angle of a sine wave to represent inspiration/expiration.
    t= linspace(1,100,50);
    f=2;
    a= cos(2*pi*f*t)

    subplot(2,4,3); cla
    x=0.01:0.0025:2;
    % ECGcomplete; % ECG specs:
    li=15/72;   % hr?
    %             li=10/50;   % hr?

    a_pwav=0.25;
    d_pwav=0.09;
    t_pwav=0.16;

    a_qwav=0.025;
    d_qwav=0.066;
    t_qwav=0.166;

    a_qrswav=1.6;
    d_qrswav=0.11;

    a_swav=0.25;
    d_swav=0.066;
    t_swav=0.09;

    a_twav=0.35;
    d_twav=0.142;
    t_twav=0.2;

    a_uwav=0.035;
    d_uwav=0.0476;
    t_uwav=0.433;
    %pwav
    pwav=p_wav(x,a_pwav,d_pwav,t_pwav,li);

    %qwav output
    qwav=q_wav(x,a_qwav,d_qwav,t_qwav,li);


    %qrswav output
    qrswav=qrs_wav(x,a_qrswav,d_qrswav,li);

    %swav output
    swav=s_wav(x,a_swav,d_swav,t_swav,li);


    %twav output
    twav=t_wav(x,a_twav,d_twav,t_twav,li);


    %uwav output
    uwav=u_wav(x,a_uwav,d_uwav,t_uwav,li);

    %ecg output
    ecg=pwav+qrswav+twav+swav+qwav+uwav;
    % to get the waveforms aligned, crop
    ecgC= imresize(ecg(25:end-120), [1,length(x)]);
    %  figure(1)
    plot(x,ecgC, 'k-','linew',2);
    legend('ECG', 'autoupdate','off','fontsize',fntsize-4)
    % xlabel('stride-cycle (%)');
    set(gca,'fontsize', fntsize, 'XTick', [0 1 2],'XTickLabel', {'0', '50','100'})
    set(gca,'ytick',[]); ylabel('a.u');
    box off
    ylim([2 5])
    hold on
    YDr= YD*6;
    plot(x, imresize(YDr, [1, length(x)]), 'b:','linew',2);
    %
    subplot(2,4,7); cla
    resp=plot(t,a, 'k-', 'linew',2);
    box off
    hold on;
    ang= plot(t,angle(hilbert(a))/2)
    shg
    ylim([-2 4])
    xlabel('stride-cycle (%)');
    set(gca,'fontsize', fntsize)
    [lg, hobh]=legend([resp,ang],{'Resp. cycle', 'Resp. phase'} ,'Location','North', 'NumColumns',2,'fontsize',fntsize-4);


    %shrink size of line obks
    linesh = findobj(hobh, 'type', 'Line');
    % restrict to 50 % size, centred on original centre.
    for ih=[1,3]
        cntr=mean(linesh(ih).XData);
        adj = diff(linesh(ih).XData)/4;
        linesh(ih).XData= [cntr-adj, cntr+adj];
    end
    set(gca,'ytick',[]); ylabel('a.u');
    set(lg,'NumColumns',2)
    box off
    cd(procdatadir)
%%
    % now plot the gabor sensitivity data:
    % load gabor data

    setmydirs_detectvGabor;

    % load data
    cd([procdatadir filesep 'GFX'])
    %load obs and shuffld data for plots:
    load('GFX_Data_inGaits.mat','GFX_TargPosData');

    % load('GFX_Data_inGaits_FourierFits.mat')
    %
    groupData = [];
    for ippant=1:size(GFX_TargPosData,1)
        % ppant, speed, ft.
        groupData(ippant,:) = GFX_TargPosData(ippant,3,3).doubgc_binned_dprime;

    end
    subplot(2,4,4);cla
    pdata= mean(groupData,1);
    plot(xvec, pdata, ['o-k']);
    hold on;

    stE= CousineauSEM(groupData);
    errorbar(xvec, pdata, stE,'k','linestyle','none')
    [f,gof]= fit(xvec',pdata',  'fourier1');
    set(gca,'color', [.9 .9 .9])
    % plot:
    hold on;
    %             yyaxis right
    Fh=plot(f, xvec, pdata);%,
    Fh(2).LineWidth = 2;
    Fh(2).Color = 'b';
    Fh(1).Visible= 'off';
    lgh=Fh(2);
    ylabel("sensitivity (d')");
    legend off
    title(['Pilot, n=' num2str(size(groupData,1))],'HorizontalAlignment','left');
    set(gca,'fontsize', fntsize, 'XTick', [0 50 100],'XTickLabel', {'0', '50','100'})
    xlabel('')

    % now the auditory task:
    set_myOnlineDirectories_AUD;
    cd(savedatadir); cd('GFX');
    load('GFX_Data_inGaits.mat','GFX_TargPosData');
    %
    groupData = [];
    for ippant=1:size(GFX_TargPosData,1)
        % ppant, speed, ft.
        groupData(ippant,:) = GFX_TargPosData(ippant,3,3).doubgc_binned_dprime;

    end
    subplot(2,4,8);cla
    pdata= mean(groupData,1);

    plot(xvec,pdata, ['o-k']);
    hold on;
    stE= CousineauSEM(groupData);
    errorbar(xvec, pdata, stE,'k','linestyle','none')
    [f,gof]= fit(xvec',pdata',  'fourier1');
    set(gca,'color', [.9 .9 .9])
    % plot:
    hold on;

    Fh=plot(f, xvec, pdata);
    Fh(2).LineWidth = 2;
    Fh(2).Color = 'm';
    Fh(1).Visible= 'off';
    lgh=Fh(2);
    ylabel("sensitivity (d')");
    legend off
    title(['Pilot, n=' num2str(size(groupData,1))],'HorizontalAlignment','left');
    set(gca,'fontsize', fntsize, 'XTick', [0 50 100],'XTickLabel', {'0', '50','100'})
    xlabel('')
    xlabel('stride-cycle (%)');

    shg
    %%
    export_fig 'DECRA2A.png'

end % job2



if job.plotFig2_ver2==1
    %% this job also needs the auditory sensitivity data loaded!

    setmydirs_detectv3;
    cd(procdatadir);
    cd(procdatadir); cd('GFX')
    load('GFX_Data_inGaits_FourierFits.mat');
    %%
    figure(1); clf; set(gcf,'color','w')
%     plot variation in accuracy and rt
    for ippant = 1:length(PFX_FourierNull);
        subplot(6,6,ippant);
        plot(PFX_FourierNull(ippant).TargetOns_doubgc_binned_rts_fitsRsq_Obs, 'r'); hold on
            plot(PFX_FourierNull(ippant).TargetOns_doubgc_binned_Acc_fitsRsq_Obs, 'b'); hold on
            title(num2str(ippant))
    end
    

    %% first plot the group result
    figure(2); clf
    set(gcf,'units','normalized','position', [.1 .1 .7 .5]);
    fntsize=12;
    useD = {targAcc*100, targRT};
    bars ={'Accuracy', 'RTs'};
    ys = {'(%)', '(sec)'};
    cols = {'b', 'r'};
    subsp= [1,2];

    for idata=1:2
        subplot(2,4,1);%subsp(idata))%[4,8])
        % accuracy?
        mdiff = round(mean(diff(pidx2)./2));
        xvec = pidx2(1:end-1) + mdiff;

        tmpD= useD{idata};
        % gM= nanmean(targAcc,1);
%remove the bars:

        % add fit:
        if idata==1
            yyaxis left
        else
            yyaxis right
        end
        pdata= mean(tmpD,1);
bh=plot(xvec, pdata, '-o', 'color', cols{idata});
hold on;
        [f,gof]= fit(xvec',pdata',  'fourier1');
hold on;
% plot black ver for legend
        Fhbl=plot(f, xvec, pdata, 'k');%,
        Fhbl(2).LineWidth = 2;
        Fhbl(2).Color = 'k';
        Fhbl(1).Visible= 'off'
%           
hold on
        Fh=plot(f, xvec, pdata);%,
        Fh(2).LineWidth = 2;
        Fh(2).Color = 'b';
        Fh(1).Visible= 'off'
%           [bh,Fh]=bar_wFourier(xvec, tmpD);
%         bh.FaceColor= [.9 .9 .9];
%         bh.EdgeColor=cols{idata}; hold on
%         bh.FaceAlpha = .2;
        gM = nanmean(tmpD);
        sdrange = max(gM) - min(gM);
        ylim([min(gM)-.5*sdrange max(gM)+.5*sdrange])
        
        legend([bh Fh(2)], {'data', 'fit'}, 'fontsize',fntsize-4);
        Fh(2).LineWidth=4;
        Fh(2).Color= cols{idata}
        
        % where to place the text?

        ylsat = get(gca, 'ylim');   
        p90 = abs(diff(ylsat))*.95 + ylsat(1);
     
        if idata==1
        
%       
        text(15, p90, 'Accuracy(%)','color','b','BackgroundColor','w','HorizontalAlignment','center','fontsize',fntsize,'fontweight','bold')
        else
%         ylabel(ys{idata}, 'units','normalized', 'Position',[1 1.175 1], 'Rotation', 0);%VerticalAlignment','top')
          text(85, p90, 'RT(sec)','HorizontalAlignment','center','Color','r','BackgroundColor','w','fontsize',fntsize,'fontweight','bold')
        end
       ylabel('')
        if idata==2
            xlabel('stride-cycle (%)');
            set(gca,'YColor', 'k')
        else
            xlabel(' ')
            set(gca,'YColor', 'k')
        end
        set(gca,'fontsize',fntsize)
        %
        hghts= [.875 , .4];
        % legend([Fh(2)], {'group fit'}, 'fontsize',fntsize-4, 'NumColumns',2, 'units','normalized', 'Position', [0.175, hghts(idata),.1 .05]);%,'Position', [0.2 0.2 .1 .1]);
        legend([Fhbl(2)], {'group fit'}, 'fontsize',fntsize-2, 'NumColumns',2, 'units','normalized', 'Location', 'South')%'Position', [0.175, hghts(idata),.1 .05]);%,'Position', [0.2 0.2 .1 .1]);

        shg

    end
    box off
    set(gcf,'color','w')
%
    %who some examples of ppant fits.
    % ppant 6 vs 29 for 2 vs 4 cps in Accuracy;
    % clf
    accpp = [5,29];

    % offset =

    % cla
    hghts= [.45, .92];
    for ipp=1:2
        subplot(2,4,2); hold on
        pdata = targAcc(accpp(ipp),:);
        % offset
        pdata= (pdata-mean(pdata)) + ipp/2;
        plot(xvec, pdata, 'ko-');
        % add fit:

        [f,gof]= fit(xvec',pdata',  'fourier1');

        % plot:
        hold on;
        %             yyaxis right
        Fh=plot(f, xvec, pdata);%,
        Fh(2).LineWidth = 2;
        Fh(2).Color = 'b';
        Fh(1).Visible= 'off'
        lgh=Fh(2);
%         text(1, ipp/2+.25, ['ppt. ' num2str(ipp)], 'fontsize', fntsize-2,'color','b');
        % save p data for next plot

        if ipp==2
            YD = Fh(2).YData;

        end
        yla= get(gca,'ylim'); % extend by 20%
        dY= yla(2) - yla(1);
        ylim([yla(1) yla(1)+dY*1.75])
        % legend(lgh,'participant fit','fontsize',fntsize-4, 'NumColumns',2, 'units','normalized', 'Position', [0.51, hghts(2),.1 .05]);%,'Position', [0.2 0.2 .1 .1]);
        legend(lgh,'participant fits','fontsize',fntsize-2, 'NumColumns',2, 'units','normalized', 'Location', 'NorthOutside');% [0.51, hghts(2),.1 .05]);%,'Position', [0.2 0.2 .1 .1]);
    end
    set(gca,'YTickLabel',[])
    ylabel('Accuracy (%)')
    set(gca,'fontsize',fntsize)
                xlabel('stride-cycle (%)');

    shg
% 
    % now plot the gabor sensitivity data:
    % load gabor data

    setmydirs_detectvGabor;

    % load data
    cd([procdatadir filesep 'GFX'])
    %load obs and shuffld data for plots:
    load('GFX_Data_inGaits.mat','GFX_TargPosData', 'pidx2');
     mdiff = round(mean(diff(pidx2)./2));
        xvec = pidx2(1:end-1) + mdiff;

    % load('GFX_Data_inGaits_FourierFits.mat')
    %
    groupData = [];
    for ippant=1:size(GFX_TargPosData,1)
        % ppant, speed, ft.
        groupData(ippant,:) = GFX_TargPosData(ippant,3,3).doubgc_binned_dprime;

    end
    subplot(2,4,3);cla
    pdata= mean(groupData,1);
    plot(xvec, pdata, ['o-k']);
    hold on;

    stE= CousineauSEM(groupData);
    errorbar(xvec, pdata, stE,'k','linestyle','none')
    [f,gof]= fit(xvec',pdata',  'fourier1');
    set(gca,'color', [.9 .9 .9])
    % plot:
    hold on;
    %             yyaxis right
    Fh=plot(f, xvec, pdata);%,
    Fh(2).LineWidth = 2;
    Fh(2).Color = 'b';
    Fh(1).Visible= 'off';
    lgh=Fh(2);
    ylabel("sensitivity (d')");
    legend off
%     title(['Pilot, n=' num2str(size(groupData,1))],'HorizontalAlignment','left');
    set(gca,'fontsize', fntsize, 'XTick', [0 50 100],'XTickLabel', {'0', '50','100'})
    xlabel('')
                xlabel('stride-cycle (%)');

    %% now the auditory task:
    set_myOnlineDirectories_AUD;
    cd(savedatadir); cd('GFX');
    load('GFX_Data_inGaits.mat','GFX_TargPosData','pidx2');
     mdiff = round(mean(diff(pidx2)./2));
        xvec = pidx2(1:end-1) + mdiff;
    %
    groupData = [];
    for ippant=1:size(GFX_TargPosData,1)
        % ppant, speed, ft.
        groupData(ippant,:) = GFX_TargPosData(ippant,3,3).doubgc_binned_dprime;

    end
    subplot(2,4,4);cla
    pdata= mean(groupData,1);

    plot(xvec,pdata, ['o-k']);
    hold on;
    stE= CousineauSEM(groupData);
    errorbar(xvec, pdata, stE,'k','linestyle','none')
    [f,gof]= fit(xvec',pdata',  'fourier1');
    set(gca,'color', [.9 .9 .9])
    % plot:
    hold on;

    Fh=plot(f, xvec, pdata);
    Fh(2).LineWidth = 2;
    Fh(2).Color = 'm';
    Fh(1).Visible= 'off';
    lgh=Fh(2);
    ylabel("sensitivity (d')");
    legend off
%     title(['Pilot, n=' num2str(size(groupData,1))],'HorizontalAlignment','left');
    set(gca,'fontsize', fntsize, 'XTick', [0 50 100],'XTickLabel', {'0', '50','100'})
    xlabel('')
    xlabel('stride-cycle (%)');

    shg
    %%
    export_fig 'DECRA2A.png'

end % job2


if job.plotFig3==1
    %%
    setdirs_subit_v1;



    cd([procdatadir filesep 'GFX']);

    load('GFX_Data_inGaits', ...
        'GFX_PropDetectData','GFX_PropDetectData_nPer',...
        'GFX_headY', 'GFX_TargPosData','GFX_RespPosData',...
        'subjIDs', 'pidx1', 'pidx2', 'gaittypes');
    %%
    % data is in the GFX_TargPosData.
    allD= [];
    for ippant = 1:size(GFX_TargPosData,1)

        %dims are speed (combined), feet( combined)
        allD(ippant,:,:) = GFX_TargPosData(ippant, 3, 3 ).doubgc_byDigitRange_binned_Acc;
    end

    figure(2);clf

    %x axis:          %approx centre point of the binns.
    mdiff = round(mean(diff(pidx2)./2));
    xvec = pidx2(1:end-1) + mdiff;
    fitcols= {'m', 'b'};
    legh=[];
    for irange = 1:2
        tmp = squeeze(allD(:,irange,:));
        plot(xvec,mean(tmp,1), 'o-k');
        hold on;
        % add fit:

        stE= CousineauSEM(tmp);
        errorbar(xvec, mean(tmp,1),stE,'k','LineStyle','none')
        [f,gof]= fit(xvec',mean(tmp,1)',  'fourier1');

        % plot:
        hold on;
        %             yyaxis right
        Fh=plot(f, xvec, mean(tmp,1));%,
        Fh(2).LineWidth = 2;
        Fh(2).Color = fitcols{irange};
        Fh(1).Visible= 'off';
        legh(irange)= Fh(2);

    end
    legend off
    ylabel('Accuracy');
    set(gca,'fontsize',fntsize)
    xlabel('stride-cycle (%)')
%     legend(legh, {'subitizing', 'estimation'}, 'NumColumns',2);
%     text(0, .96, 'Pilot data, \itn\rm=17', 'fontsize', fntsize)
    set(gcf,'color','w')

    text(0, .94, 'Subitizing', 'fontsize',fntsize, 'color','m','fontweight','bold')

    text(0, .86, 'Estimation', 'fontsize',fntsize, 'color','b','fontweight','bold')

    % shrink legend components
      % [lg,hh]=legend(legh, {'\bf\theta band', '\bf\alpha band'})
%           linesh = findobj(hh, 'type', 'Line');
    % restrict to 50 % size, centred on original centre.
%     for ih=[1,3]
%         cntr=mean(linesh(ih).XData);
%         adj = diff(linesh(ih).XData)/10;
%         linesh(ih).XData= [cntr-adj, cntr+adj];
%     end
end
%% -------------------------------------
if job.plotFig4_EEG_mobi==1
    %%
    % GRANT_plotEEGsummary
    % here we will make a fancy EEG figure to promote further experiments.

    % load EEG data:
    % cd('C:\Users\mdav0285\mobiThings\Processed');
    % cd('C:\Users\mdav0285\mobiThings\Processed');
    cd('C:\Users\mdav0285\Documents\Data\mobiThings\Processed')
    procdatadir=pwd;
    % load data (EEG)
    cd(procdatadir)
    cd('GFX');
    load('GFX_gaitEpochedData_summary');

    %%
    clf
    DSI_chans={'Fp1', 'Fp2',...
        'F7','F3','Fz', 'F4', 'F8', ...
        'T3','C3', 'Cz', 'C4','T4',...
        'T5','P3', 'Pz','P4','T6','O1', 'O2'};

    speedCols={'b',[1, 171/255, 64/255], 'k'}; % blue, "normal" yellow

%     frontChans = {'Fp1', 'Fp2', 'F3', 'Fz', 'F4', 'C3', 'Cz', 'C4'};
        frontChans = { 'F3', 'Fz', 'F4', 'C3', 'Cz', 'C4'};
    frontChan_IDX = find(contains(DSI_chans, frontChans));

    backChans = {'P3', 'Pz', 'P4', 'O1', 'O2', 'T5', 'T6'};
    backChan_IDX = find(contains(DSI_chans, backChans));


    % also load Head data for overlay:
    setmydirs_detectv3
    % next plot overall accuracy and RT,
    cd([procdatadir filesep 'GFX']);
    load('GFX_Data_inGaits.mat', 'GFX_headY');

    fntsize= 12;
    %concat:
    GFX_plotHead=[];

    for ippant = 1:length(GFX_headY);

        GFX_plotHead(ippant,:) = GFX_headY(ippant).gc;
        pidx=ceil(linspace(1,100,11)); % length n-1

    end

    % Prepare data for Time frequency:
    %prepare subj level data:


    % TF within channels. (within channels and both speeds).
    frontTF_chanAv = squeeze(nanmean(nanmean(GFX_TF_EEG(:,:,frontChan_IDX,:,:),3),2));
    backTF_chanAv=squeeze(nanmean(nanmean(GFX_TF_EEG(:,:,backChan_IDX,:,:),3),2)); %

    %or all channels by speed:
    % norm of TF over all channels.
    slowTF_chanAv = squeeze(nanmean(GFX_TF_EEG(:,1,:,:,:),3));
    naturalTF_chanAv=squeeze(nanmean(GFX_TF_EEG(:,2,:,:,:),3));

    diff_TFspeed = naturalTF_chanAv- slowTF_chanAv;
% diff_TFspeed = slowTF_chanAv;
        
diff_TFloc = frontTF_chanAv-backTF_chanAv;


    plotD= {diff_TFspeed, diff_TFloc};

    % set up plot
    figure(1);clf;
    set(gcf,'color', 'w', 'units', 'normalized', 'position', [.1 .1 .6 .4]);
    shleg=[];

    stepAxis = [1:260] ./2;

    fAt = [1,2,5,10,15,30];
    titlesAre= {'Normal - Slow walk', 'Frontal - Occipital'};
    % titlesAre= {'Normal - Slow walk'};

    %

    clf
    % different plot levels to show sig:

    plotLevels= {[-4, 2.15],[-3.6, 2.15]};
figure(1); clf; set(gcf,'color','w')
    for id= 1%:2
        figure(1)
subplot(3,2,id);

subplot(2,3,id);  % for colorbar
        hold on;
%         % note in this version we need to pad the head.
%         headH=mean(GFX_plotHead,1);
%         headPad = [nan(1,30), imresize(headH,[1,200]), nan(1,30)];
%         lp=plot(stepAxis, headPad, 'mo-', 'LineWidth', 2,'MarkerIndices',1:5:260);
%         set(gca,'ytick', [], 'YColor', 'w');
%         xlim([stepAxis(1) stepAxis(end)])
% %         ylabel('Head height (m)')
%        cfg=[];
%        cfg.gap=0;
%        cfg.heightMod=2;
%        subplotBelow(gca, cfg)

%
        diffTF= plotD{id};

        diffM = squeeze(nanmean(diffTF,1));
        %     title({'Frontal vs Occipital EEG channels'}); hold on;
        %
        contourf(stepAxis, frex, diffM*100, 10, 'LineStyle',':');
        hold on;
        contour(stepAxis, frex, diffM*100, plotLevels{id}, 'k-', 'LineWidth', 4)

        %     imagesc(stepAxis, frex, diffM);
        xlabel('step-cycle (%)');

        if id==1
            c=colorbar;
            ylabel(c, '% power change', 'fontsize', fntsize);
%             colormap('viridis')
rdbu = cbrewer('div', 'RdYlBu',100);
rdbu(rdbu>1)=1;
rdbu(rdbu<0)=0;
colormap(flipud(rdbu));
         ylabel('Frequency (Hz)')
        else
            ylabel('Frequency (Hz)')
        end
        set(gca,'fontsize',fntsize, 'YScale', 'log', 'ytick', [fAt],...
            'YTickLabel', {'1', '2', '5', '10', '15',  '30', '40'}, ...
            'Xtick', [15,65, 115], 'XTickLabel', {'1', '50', '100'});
        if id==2
%             set(gca,'YTickLabel', []);
        end
        caxis([-20 20]);
%         caxis([-8 8])
%         title(titlesAre{id})
%         axis square
        %
        %stats
        pvals=zeros(size(diffTF,2), size(diffTF,3));
        for ifreq= 1:size(diffTF,2)
            for isamp= 1:size(diffTF,3)

                [H,pval] = ttest(diffTF(:,ifreq,isamp));
                if pval<.05
                    pvals(ifreq, isamp) = pval;
                end
            end
        end

        %%
figure(1);
subplot(2,3,2);cla
yyaxis right
 headH=mean(GFX_plotHead,1);
        headPad = [nan(1,30), imresize(headH,[1,200]), nan(1,30)];
        lp=plot(stepAxis, headPad, 'm-', 'LineWidth', 2);
        lp.Color = [1,0,1,.2];
        set(gca,'ytick', [], 'YColor', 'k');
        xlim([stepAxis(1) stepAxis(end)])

yyaxis left
% plot activity in bands:
bands = [1,5; 8,12];
powcols= {'b','r'};
legh=[];
for iband = 1:2
    % tmpP
    frange = dsearchn(frex', bands(iband,:)');
    tmpPower = squeeze(nanmean(diffTF(:, frange(1):frange(2),:),2))*100;
    pM = mean(tmpPower,1);
    pSte= CousineauSEM(tmpPower);

%     plot(stepAxis, pM, powcols{iband}); hold on;
    sh=shadedErrorBar(stepAxis, pM, pSte,{powcols{iband}, 'linew',2},1); hold on
legh(iband)=sh.mainLine;
end
text(85, 35, ['\bf\theta (2-7 Hz)'],'color','r','fontsize',10); 
text(85, 29, ['\bf\alpha (8-12 Hz)'],'color','b','fontsize',10); 

        xlabel('step-cycle (%)');
          %shrink size of line obks
    % [lg,hh]=legend(legh, {'\bf\theta band', '\bf\alpha band'})
%           linesh = findobj(hh, 'type', 'Line');
    % restrict to 50 % size, centred on original centre.
%     for ih=[1,3]
%         cntr=mean(linesh(ih).XData);
%         adj = diff(linesh(ih).XData)/10;
%         linesh(ih).XData= [cntr-adj, cntr+adj];
%     end
    ylim([-30 40])
xlim([stepAxis(1) stepAxis(end)]);
%
 set(gca,'fontsize',fntsize, 'Xtick', [15,65, 115], 'YColor', 'k','XTickLabel', {'1', '50', '100'});
%%

        figure(3);
        subplot(2,3,id)
       
            topoplot([], elocs, 'emarker', {'o','k',5,5})
            hold on;
            topoplot(zeros(1,19), elocs, 'emarker2', {frontChan_IDX, 'o','r'})
          topoplot(zeros(1,19), elocs, 'emarker2', {backChan_IDX, 'o','b'})
          
            set(gcf,'color','w')
subplot(2,3,id+2); % plot the colour bar we need
%             contour(stepAxis, frex, diffM*100, plotLevels{id}, 'w-', 'LineWidth', 2)

        % change topo to white.
        cmap = colormap;
        cmap= [1,1,1];
        colormap(cmap)
%         figure(1); hold on
    end % id
end
%%
if job.plotACNSpanel1== 1; % the detectv3 summary, new spacing



    %% this job also needs the auditory sensitivity data loaded!

    setmydirs_detectv3;
    cd(procdatadir);
    cd(procdatadir); cd('GFX')
    load('GFX_Data_inGaits_FourierFits.mat');
    %%
    figure(1); clf; set(gcf,'color','w')
%     plot variation in accuracy and rt
    for ippant = 1:length(PFX_FourierNull);
        subplot(6,6,ippant);
        plot(PFX_FourierNull(ippant).TargetOns_doubgc_binned_rts_fitsRsq_Obs, 'r'); hold on
            plot(PFX_FourierNull(ippant).TargetOns_doubgc_binned_Acc_fitsRsq_Obs, 'b'); hold on
            title(num2str(ippant))
    end
    

    %% first plot head height change, then group result
     figure(2); clf

    set(gcf,'units','normalized','position', [0 0 .4 .4]);
    subplot(2,2,[1,3]); cla

    % yyaxis right;
    % tmp = nanmean(targCounts,1);
    bh=bar(1:40, nanmean(targCounts,1));
    bh.FaceColor= [.9 .9 .9];
    bh.EdgeColor=[.5 .5 .5]; hold on
    bh.FaceAlpha = .2;
    ylim([0 20])
    set(gca,'YColor', 'k','fontsize',12);
    ylabel('Target counts');%, 'HorizontalAlignment','right')
    % yyaxis left
    plot([1:40], imresize(mean(headY,1),[1,40])*100+20, ['ko-'],'MarkerSize',10,'MarkerIndices',1:1:200);
    ylim([-4 25])
    hold on
    ht= plot(nan, nan, ['ko'],'MarkerSize',10,'MarkerIndices',1:2:200)
    xlabel('');% stride-cycle');
    legend(ht, 'Head height','fontsize', fntsize, 'units','normalized','Position', [0.196 .89 .21 .07])
    % ylabel(['Head height (m)'],'fontsize',13,'HorizontalAlignment','left');
    set(gca,'ytick',[0:2:18],'XTickLabel',[],'YColor','k', 'fontsize', fntsize);
    box on
    %
    %
    %
    hold on
    % subplot(2,2,3);
    % plot([1:200]./2, repmat(1,[1,200]), 'w')
    xlabel('stride-cycle (%)');
    % set(gca,'ytick',[], 'YColor', 'w')

    
    fntsize=12;
    useD = {targAcc*100, targRT};
    bars ={'Accuracy', 'RTs'};
    ys = {'(%)', '(sec)'};
    cols = {'b', 'r'};
    subsp= [1,2];

    for idata=1:2
        subplot(2,2,[2,4]);%subsp(idata))%[4,8])
        % accuracy?
        mdiff = round(mean(diff(pidx2)./2));
        xvec = pidx2(1:end-1) + mdiff;

        tmpD= useD{idata};
        % gM= nanmean(targAcc,1);
%remove the bars:

        % add fit:
        if idata==1
            yyaxis left
        else
            yyaxis right
        end
        pdata= mean(tmpD,1);
bh=plot(xvec, pdata, '-o', 'color', cols{idata});
hold on;
        [f,gof]= fit(xvec',pdata',  'fourier1');
hold on;
% plot black ver for legend
        Fhbl=plot(f, xvec, pdata, 'k');%,
        Fhbl(2).LineWidth = 2;
        Fhbl(2).Color = 'k';
        Fhbl(1).Visible= 'off'
%           
hold on
        Fh=plot(f, xvec, pdata);%,
        Fh(2).LineWidth = 2;
        Fh(2).Color = 'b';
        Fh(1).Visible= 'off'
%           [bh,Fh]=bar_wFourier(xvec, tmpD);
%         bh.FaceColor= [.9 .9 .9];
%         bh.EdgeColor=cols{idata}; hold on
%         bh.FaceAlpha = .2;
        gM = nanmean(tmpD);
        sdrange = max(gM) - min(gM);
        ylim([min(gM)-.25*sdrange max(gM)+.5*sdrange])
        ylsat = get(gca, 'ylim');
        legend([bh Fh(2)], {'data', 'fit'}, 'fontsize',fntsize-4);
        Fh(2).LineWidth=4;
        Fh(2).Color= cols{idata}
%         if idata==1
        
        ylabel(ys{idata})%, 'units','normalized', 'Position',[-.1 1 1], 'Rotation', 0);%VerticalAlignment','top')
        
            
%         ya= get(gca,'ylim');
        % place text at 95% mark.
%         yAt = (ya(2) - ya(1))*1.1 + ya(1);
%         text(5, yAt, bars{idata}, 'fontsize',fntsize)
        % text(bars{idata}, 'HorizontalAlignment','right')


        if idata==2
            xlabel('stride-cycle (%)');
            set(gca,'YColor', 'k', 'YTick', .35:.01:.4)
            % move RT label to fit:
%         ylabel([bars{idata} ' '  ys{idata} ], 'units','normalized', 'Position',[1.1 .6 1], 'FontWeight', 'bold','Color','r')
        else
            xlabel(' ')
            set(gca,'YColor', 'k')
%                        ylabel([bars{idata} ' '  ys{idata} ], 'units','normalized', 'Position',[-.15 .6 1], 'FontWeight', 'bold','Color','b')
        end
        set(gca,'fontsize',fntsize)
        %
        hghts= [.875 , .4];
        % legend([Fh(2)], {'group fit'}, 'fontsize',fntsize-4, 'NumColumns',2, 'units','normalized', 'Position', [0.175, hghts(idata),.1 .05]);%,'Position', [0.2 0.2 .1 .1]);
        legend([Fhbl(2)], {'group fit'}, 'fontsize',fntsize, 'NumColumns',2, 'units','normalized', 'Location', 'NorthEast')%'Position', [0.175, hghts(idata),.1 .05]);%,'Position', [0.2 0.2 .1 .1]);

        shg
        textheight = ylsat(1) + (ylsat(2)-ylsat(1))*.8;
        text(1, textheight ,'Accuracy (%)','Fontweight','bold','color','b','fontsize',fntsize)

        text(xvec(end), textheight ,'RT (sec)','Fontweight','bold','color','r','fontsize',fntsize,'HorizontalAlignment','right')

    end
    set(gca,'xtick',[])
    box on
    set(gcf,'color','w')
end
