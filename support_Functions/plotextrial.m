
function plotextrial(cfg)

setmydirs_detectv3;
cd(procdatadir);
subjID= 'ABB';
itrial= 124;
cd(procdatadir)
lfile  = dir([pwd filesep subjID '*']);
load(lfile.name)
fntsize= cfg.fntsize;

%
plot(trialInfo(124).times, HeadPos(124).Y, 'k', 'LineWidth',2);
% axis tight;
xlabel('Trial time (sec)');
ylabel('Head height (m)');
set(gca,'fontsize', fntsize);
ylim([1.69 1.78])
%
% tOnsets = HeadPos(124).tr   
relT= find(trial_summaryTable.trial == itrial);
tOnsets = trial_summaryTable.targOnset(relT);
tResult = trial_summaryTable.targCor(relT);
hold on;
trialTimes = trialInfo(itrial).times;
 HeadData= HeadPos(itrial).Y;
txtHeight = 1.7;
 for itarg = 1:length(tOnsets)
   
     pAt = dsearchn(trialTimes, tOnsets(itarg)');

if tResult(itarg)==1% Hit.
%     ht=text(tOnsets(itarg), txtHeight, 'H', 'HorizontalAlignment', 'center')
    hl=plot(tOnsets(itarg), HeadData(pAt), 'ro', 'MarkerFaceColor', 'r');
    
    plot(tOnsets(itarg), 1.695, 'vk', 'MarkerFaceColor', 'k') 

else % Miss

%     mt=text(tOnsets(itarg), txtHeight, 'M', 'HorizontalAlignment','center')
    ml=plot(tOnsets(itarg), HeadData(pAt), 'ro', 'MarkerFaceColor', 'w');
    plot(tOnsets(itarg), 1.695, 'vk', 'MarkerFaceColor', 'w') 

end

% connect to base
plot([tOnsets(itarg) tOnsets(itarg)], [1.695 1.72], 'k-')
% plot(tOnsets(itarg), 1.7, 'dk', 'MarkerFaceColor', 'k') 
 end
 legend([hl,ml], {'Hit', 'Miss'})

 text(0.1, 1.71, {['Target'];['onset:']}, 'HorizontalAlignment', 'left', 'fontsize', fntsize-2);
 xlim([0 9])
 box off

end

