%GRANT_plotProbeSummary_v2

% This script will produce a document ready figure of our probe paradigm.

% 1 panel showing an example trial of head positiond ata (Y axis), with target onsets.
% 2 panel showing epoched version with target accuracy on Y.
% 3 rd panel showing RT per location.


%This version, omits the full trial head-profile, using a stride-cycle
%instead (to show step portions in the grant).



% example trial( good one, is participant ABB, trial 124);


% load Head Y data. 


setmydirs_detectv3;
cd(procdatadir)
% show ppant numbers:
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)

%
fntsize = 12;


subjID= 'ABB';
itrial= 124;
lfile  = dir([pwd filesep subjID '*']);
load(lfile.name)


%% 
figure(1);
% clf; set(gcf, 'color', 'w', 'units', 'normalized', 'position', [.1 .1  .3 .6]);
 set(gcf, 'color', 'w', 'units', 'normalized', 'position', [.1 .1  .65 .25]);
%%
 clf;

nGaits=2;
addFit=1;
subplot(1,3,1);
% subplot(2,2,1:2);


% %
% plot(trialInfo(124).times, HeadPos(124).Y, 'k');
% % axis tight;
% xlabel('Trial time (sec)');
% ylabel('Head height (m)');
% set(gca,'fontsize', fntsize);
% ylim([1.69 1.78])
% %
% % tOnsets = HeadPos(124).tr   
% relT= find(trial_summaryTable.trial == itrial);
% tOnsets = trial_summaryTable.targOnset(relT);
% tResult = trial_summaryTable.targCor(relT);
% hold on;
% trialTimes = trialInfo(itrial).times;
%  HeadData= HeadPos(itrial).Y;
% txtHeight = 1.7;
%  for itarg = 1:length(tOnsets)
%    
%      pAt = dsearchn(trialTimes, tOnsets(itarg)');
% 
% if tResult(itarg)==1% Hit.
% %     ht=text(tOnsets(itarg), txtHeight, 'H', 'HorizontalAlignment', 'center')
%     hl=plot(tOnsets(itarg), HeadData(pAt), 'ro', 'MarkerFaceColor', 'r');
%     
%     plot(tOnsets(itarg), 1.695, 'vk', 'MarkerFaceColor', 'k') 
% 
% else % Miss
% 
% %     mt=text(tOnsets(itarg), txtHeight, 'M', 'HorizontalAlignment','center')
%     ml=plot(tOnsets(itarg), HeadData(pAt), 'ro', 'MarkerFaceColor', 'w');
%     plot(tOnsets(itarg), 1.695, 'vk', 'MarkerFaceColor', 'w') 
% 
% end
% 
% % connect to base
% plot([tOnsets(itarg) tOnsets(itarg)], [1.695 1.72], 'k-')
% % plot(tOnsets(itarg), 1.7, 'dk', 'MarkerFaceColor', 'k') 
%  end
%  legend([hl,ml], {'Hit', 'Miss'})
% 
%  text(0.1, 1.71, {['Target'];['onset:']}, 'HorizontalAlignment', 'left', 'fontsize', fntsize-2);
%  xlim([0 9])
%  box off
% tightfig;


% next plot overall accuracy and RT,
  cd([procdatadir filesep 'GFX']);
    load('GFX_Data_inGaits.mat');

%concat:
GFX_plotHead=[];
GFX_gAcc=[];
GFX_grts=[];
for ippant = 1:size(GFX_TargPosData,1);


    if nGaits==1
GFX_plotHead(ippant,:) = GFX_headY(ippant).gc;
GFX_gAcc(ippant,:)=GFX_TargPosData(ippant,3).gc_binned_Acc;
GFX_grts(ippant,:)=GFX_TargPosData(ippant,3).gc_binned_rts;
pidx=ceil(linspace(1,100,21)); % length n-1
    else
        GFX_plotHead(ippant,:) = GFX_headY(ippant).doubgc;
GFX_gAcc(ippant,:)=GFX_TargPosData(ippant,3).doubgc_binned_Acc;
GFX_grts(ippant,:)=GFX_TargPosData(ippant,3).doubgc_binned_rts;
pidx=ceil(linspace(1,200,41)); % length n-1

    end

end



%figspecs:

mdiff = round(mean(diff(pidx)./5));
xvec = pidx(1:end-1) + mdiff;
useD = {GFX_gAcc.*100, GFX_grts};
usecols= {'b', 'r'};
ylabs = {'Accuracy (%)', 'RT (sec)'};

if nGaits==1
yticksat = {64:2:80, .34:.005:38};
else
    yticksat = {64:2:80, .34:.01:38};
end

% plot head data first
% % % % % % Head DATA:
clf
subplot(3,3,[1, 4])
plotHead = GFX_plotHead;
% plotHeadn = plotHead./ mean(gM);
%             rescale plotHead, so that the max, is above the ylim.
plotHeadM = mean(plotHead,1 ,'omitnan');
%             norm to 1.
% plotHeadM = plotHeadM./max(plotHeadM);
%adjust so max is (just) above ylim.
% ylimsn = get(gca, 'YLim');
% plotHeadM= plotHeadM.* (ylimsn(2) + .1*diff(ylimsn));

% yyaxis right
stEH= CousineauSEM(plotHead);
pH=shadedErrorBar(1:size(plotHead,2),plotHeadM, stEH,'k',1)

% yrange = max(plotHeadM) - min(plotHeadM);
%             ylim([(min(plotHeadM)- 1.2*yrange), max(plotHeadM)+.1*yrange]);

% if nGaits==1
%     ylim([(min(plotHeadM)) max(plotHeadM)+.2*yrange]);
% else
%     ylim([(min(plotHeadM))-1.7*yrange, max(plotHeadM)+.3*yrange]);
% 
% end
hold on
ylabel('Head height [m]');
hold on
pH=plot(xvec, plotHeadM(xvec), 'ok' )

set(gca,'xtick', [], 'ytick', [],'YColor', 'k','fontsize', fntsize);
% and a blank for the x axis:
subplot(3,3,7)
plot(xvec,nan(1,length(xvec)));
box off
% AXES GENERAL
axis tight
set(gca,'ytick', [])
midp=xvec(ceil(length(xvec)/2));
xlim([0 xvec(end)])
set(gca,'fontsize', fntsize, 'xtick', [1, midp, xvec(end)], 'XTickLabels', {'0', '50', '100'})
xlabel('% gait-cycle');
shg
%

shg
%%
pspots= [2,8; 3,9];
underspots= [8, 9]; % for underneath
for idata= 1:2
% subplot(1,3,1+idata)
% subplot(2,2,2+idata)
subplot(3,3,pspots(idata,:))


%% % % % % % BAR DATA:
ppantData = useD{idata};
gM = squeeze(mean(ppantData, 'omitnan'));
stE = CousineauSEM(ppantData);
% finely sampled bar, each gait "%" point.
bh=bar(xvec, gM);
bh.EdgeColor= usecols{idata}
hold on;
errorbar(xvec, gM, stE, ...
    'color', 'k',...
    'linestyle', 'none',...
    'linew', 1);
bh.FaceColor = usecols{idata};

% clip y labels.
set(gca,'Ytick', yticksat{idata})
bh.FaceAlpha= .2;
legp(idata)= bh;

%adjust view to capture most of the variance:
sdrange = max(gM) - min(gM);
if nGaits==1
ylim([min(gM)-.5*sdrange max(gM)+.7*sdrange])
else
ylim([min(gM)-.4*sdrange max(gM)+.9*sdrange])
end
% or symmetrical: 
ylim([min(gM)-.4*sdrange max(gM)+.4*sdrange])


  %  add best fourier fit
if addFit==1           
            %             best fit (use model params below):
            [f,gof]= fit(xvec',gM',  'fourier1'); %unbounded
            
            % plot:
            hold on;
            h=plot(f, xvec, gM);%,
            h(2).LineWidth = 2;
            h(2).Color = 'k';
            % include period and Rsquared
            %treat max xvec as our full 'period'
            Hzapp = xvec(end)/ (2*pi/(f.w));
            %
%             legdetails = [sprintf('%.2f', Hzapp) ' Hz_G_C, R^2 = ' sprintf('%.2f', gof.rsquare) ];
%             legend(h(2), legdetails, 'fontsize', 15, 'autoupdate', 'off')

            h(2).Color = usecols{idata};

end
legend off
ylabel(ylabs{idata})

% % % % % % Head DATA:

% plotHead = GFX_plotHead;
% plotHeadn = plotHead./ mean(gM);
% %             rescale plotHead, so that the max, is above the ylim.
% plotHeadM = mean(plotHead,1 ,'omitnan');
% %             norm to 1.
% plotHeadM = plotHeadM./max(plotHeadM);
% %adjust so max is (just) above ylim.
% ylimsn = get(gca, 'YLim');
% plotHeadM= plotHeadM.* (ylimsn(2) + .1*diff(ylimsn));
% 
% yyaxis right
% stEH= CousineauSEM(plotHead);
% pH=shadedErrorBar(1:size(plotHead,2),plotHeadM, stEH,'k',1)
% 
% yrange = max(plotHeadM) - min(plotHeadM);
% %             ylim([(min(plotHeadM)- 1.2*yrange), max(plotHeadM)+.1*yrange]);
% 
% if nGaits==1
%     ylim([(min(plotHeadM)) max(plotHeadM)+.2*yrange]);
% else
%     ylim([(min(plotHeadM))-1.7*yrange, max(plotHeadM)+.3*yrange]);
% 
% end
% % add circles?
% hold on
% pH=plot(xvec, plotHeadM(xvec), 'ok' )
% set(gca,'ytick', [],'YColor', 'w');
% subplot(3,3,underspots(idata))
plot(xvec,nan(1,length(xvec)));
box off
% AXES GENERAL
set(gca, 'xtick',[])
midp=xvec(ceil(length(xvec)/2));
set(gca,'fontsize', fntsize, 'xtick', [1, midp, xvec(end)], 'XTickLabels', {'0', '50', '100'})
xlabel('Target onset (% gait)');
shg
% %
% 
% % AXES GENERAL
% xlabel('Target onset (% gait)');
% midp=xvec(ceil(length(xvec)/2));
% set(gca,'fontsize', fntsize, 'xtick', [1, midp, xvec(end)], 'XTickLabels', {'0', '50', '100'})
% legend off
% % lg=legend(pH, 'Head position', 'Orientation','horizontal','Location','north');
% % lg.Position = [lg.Position(1), lg.Position(2)+.05, lg.Position(3), lg.Position(4)]
% box off
end
%%
shg
% tightfig