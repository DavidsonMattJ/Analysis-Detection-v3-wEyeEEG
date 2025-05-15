% plot_exampleGaitVariability

%%  show the gait variability.
% load the timeseries of all single (or double) steps. Overlay an example
% participant (all steps), and the average for all participants.

setmydirs_detectv3
cd(procdatadir);

pfols = dir([pwd filesep '*_summary_data.mat']);
job=[];
job.plotSuppFigure=0;
job.plotReviewfigure=0;

job.plotFig1_subpanels=1; % request from Rvwr 2, to include variability in Fig1.

%% good ppant exemplars
% 5,6,7, 19,20,25, 27, 29
if job.plotSuppFigure==1
figure(1); clf;
set(gcf,'color', 'w', 'units', 'normalized', 'position', [.1 .1 .8 .6])
cd(procdatadir);
ippant = 7;
% ippant = 19;
% ippant=29
    load(pfols(ippant).name, 'doubgait_ts_raw');
    doubgait_ts_raw(doubgait_ts_raw==0) =nan;
    
   figure(1); clf;
   subplot(3,3,[1,4])
   plot(doubgait_ts_raw', 'color', [0 0 0 .1]);
%    % plot all end points?
allEnds=[];
for istep = 1:size(doubgait_ts_raw,1)
    
    endAt= find(~isnan(doubgait_ts_raw(istep,:)), 1,'last');
    hold on
%     plot(endAt, doubgait_ts_raw(istep, endAt), 'o', 'color', [1,0,0,.1]);
allEnds(istep)=endAt;
    
end
title('Example participant');
set(gca,'fontsize', 12)
ylim([-.06 .06])

ylabel('detrended head height (m)');
%change x axis.
currxtic = get(gca,'xtick');
set(gca,'xtick', 0:18:119,'xticklabel', [0:18:119]./90);
xlabel('Time (sec)')
%% 
subplot(3,3,7)
histogram(allEnds./90, 40); 
ylabel('Stride count')
xlabel('Stride duration (sec)')
text(.85,150*.9, ['(n=' num2str(length(allEnds)) ')'], 'fontsize',12)
set(gca,'fontsize', 12)
xlim([.8 1.5]);
ylim([0 150]);
box off
%%

% create supersubject
 subplot(3,3,[2 5]); cla
 %retain all Ends:
 endAt= [];
for ippant= 1:length(pfols);
    cd(procdatadir);
    load(pfols(ippant).name, 'doubgait_ts_raw');
    doubgait_ts_raw(doubgait_ts_raw==0) =nan;
    
   
   plot(doubgait_ts_raw', 'color', [0 0 0 .1]);
    hold on
% % plot all end points?
for istep = 1:size(doubgait_ts_raw,1)
    
    tmpendAt= find(~isnan(doubgait_ts_raw(istep,:)), 1,'last');
    
    endAt= [endAt, tmpendAt];
    
end
end
%%
title('All participants (supersubject)');
ylabel('detrended head height (m)');
%%
%change x axis.
currxtic = get(gca,'xtick');
set(gca,'xtick', 0:18:119,'xticklabel', [0:18:119]./90,'fontsize',12);
xlabel('Time (sec)')
xlim([0 119]);
%%
ylim([-.06 .06])
%%
subplot(3,3,8)
histogram(endAt./90); 
ylabel('Stride count')
xlabel('Stride duration (sec)')
text(.85,2100*.9, ['(n=' num2str(length(endAt)) ')'], 'fontsize',12)
set(gca,'fontsize', 12)
ylim([0 2100]);
% xlim([.8 1.6])
xlim([.8 1.5])
box off
%%
cd(procdatadir)
cd('GFX');

load('GFX_Data_inGaits', 'GFX_headY');
 % show all sub averages.
 
 subplot(3,3,[3 6]); cla
 allMin=[]
 for ippant = 1:size(GFX_headY,2);
     
     tmp =GFX_headY(ippant).doubgc_raw;
     tmp(tmp==0)=nan;
     
     % end plot at minimum point.
     [minIs, minAt] = min(tmp(80:end));
     minAt= minAt+79;
     plot(tmp(1:minAt), 'k');
     hold on
     allMin(ippant)=minAt;
 end
 ylim([-.06 .06]);
 title('Average per participant')
 ylabel('detrended head height (m)');

 set(gca,'xtick', 0:18:119,'xticklabel', [0:18:119]./90);
xlabel('Time (sec)')
set(gca,'fontsize', 12)
%%
 subplot(3,3,9)
histogram(allMin./90, 20); 
ylim([0 8])
text(.85,8*.9, ['(\itN\rm=' num2str(length(allMin)) ')'], 'fontsize',12)
ylabel('Participant count')
xlabel('Stride duration (sec)')
set(gca,'fontsize', 12)
xlim([.8 1.5])
box off
 
end

if job.plotReviewfigure==1

%% plot example of dramatic restretching:
clf;
%
subplot(321);

ippant=4;
tmp =GFX_headY(ippant).doubgc_raw;
tmp(tmp==0)=nan;
% end plot at minimum point.
     [minIs, minAt] = min(tmp(80:end));
     minAt= minAt+79;
     newtmp=tmp(1:minAt);
     plot([1:length(newtmp)]./90,newtmp, 'k:','linew', 2);
     p1l= minAt;
     hold on;
 ippant= 10    ;
 tmp =GFX_headY(ippant).doubgc_raw;
tmp(tmp==0)=nan;
% end plot at minimum point.
     [minIs, minAt] = min(tmp(80:end));
     minAt= minAt+79;
     
     newtmp = imresize(tmp(1:minAt), [1,150]);
     p2l= 150;
     plot([1:length(newtmp)]./90,newtmp, 'k-','linew', 2);
     ylabel('detrended head height');
% %      legend('participant 1', 'participant 2')
%      xlabel('stride duration (sec)');
     set(gca,'fontsize', 15)
     title('Raw')
     % plot both again (resampled).
     subplot(3,2,2); cla
     ls={':','-'};
     ic=1
     
     for ippant = [4,10]
         tmp =GFX_headY(ippant).doubgc;
         tmp= imresize(tmp,[1,100]);
         plot(tmp, ['k' ls{ic}], 'linew',2);
         hold on
         ic=ic+1;
     end
     title('stride resampled');
     set(gca,'fontsize', 15,'xtick', [0 25 50 75 100], 'xticklabel', {'0', '25%', '50%', '75%', '100%'});
%      xlabel('% stride-cycle completion');
     
     % 
     subplot(3,2,3)
     % add a hypothetical sinewave underneath.
     amps= [1, .75];
t=linspace(0,p1l,p1l)
f=2;
Y1= amps(1)*cos(2*pi*f*t)
     plot(t./90,Y1, 'r:', 'linew', 2)
     hold on;
     %
     t=linspace(0,p2l,p2l)
f=2;
Y2= amps(2)*cos(2*pi*f*t)
     plot(t./90,Y2, 'r-', 'linew',2)
     hold on;
     legend('f1', 'f2');
     ylabel('change in RT')
%      xlabel('stride duration (sec)');
set(gca,'fontsize', 15)

title('raw oscillations: f1 > f2')
%
     %\% now resample both.
       subplot(3,2,4); cla
     % add a hypothetical sinewave underneath.
    Y1re= imresize(Y1,[1,100])
    Y2re= imresize(Y2,[1,100])
    
     plot(1:100,Y1re, 'r:', 'linew', 2)
     hold on;
     plot(1:100,Y2re, 'r-', 'linew',2)
%      Y2=circshift(Y,5);
     hold on;
    av= plot(1:100, mean([Y1re;Y2re]), 'b', 'linew',2);
     legend(av, 'average');
     ylabel('change in RT')
%      xlabel('stride duration (sec)');
set(gca,'fontsize', 15)
xlim([0 100])
title('f1>f2 resampled')
set(gca,'fontsize', 15,'xtick', [0 25 50 75 100], 'xticklabel', {'0', '25%', '50%', '75%', '100%'});
%      xlabel('% stride-cycle completion');
% % now an alternative where we start with the clock phase, and resample
 

     subplot(3,2,5); cla
     % add a hypothetical sinewave underneath, now matched frequencies

     %
     t=linspace(0,p2l,p2l)
     % first just plot this sine wave up until the end of first stride:
     
f=2;
Y= amps(1)*cos(2*pi*f*t)
     plot(t(1:p1l)./90,Y(1:p1l), ':','color', [1,0,0,.5] ,'linew',2)
    
     Y1= amps(1)*Y(1:p1l);
     hold on;
     % now finish with other stride:
      plot([1:150]./90,Y2, '-','color', [1,0,0,.5], 'linew',2)
%       Y2= circshift(Y,2);
      %
     legend('f1', 'f2');
     ylabel('change in RT')
     xlabel('stride duration (sec)');
set(gca,'fontsize', 15)
title('raw if f1=f2')
%%

     %% now resample both.
       subplot(3,2,6); cla
     % add a hypothetical sinewave underneath.
     plot(1:100,imresize(Y1,[1,100]), ':','color', [1,0,0,.5], 'linew',2)
     y1re=imresize(Y1,[1,100]);
     y2re=imresize(Y2,[1,100]);
     %
     hold on;
     plot(1:100,imresize(Y2,[1,100]), '-','color', [1,0,0,.5], 'linew',2)
     
     av=plot(1:100, mean([y1re;y2re]),'b', 'linew',2);
     %%
     title('f1=f2 resampled')
     legend(av, 'Average');
     ylabel('change in RT')
     xlabel('stride duration (sec)');
set(gca,'fontsize', 15)
xlim([0 100])
% title('performance oscillations: f1 = f2')
set(gca,'fontsize', 15,'xtick', [0 25 50 75 100], 'xticklabel', {'0', '25%', '50%', '75%', '100%'});
     xlabel('% stride-cycle completion');
 
 
 
end
if job.plotFig1_subpanels==1; % request from Rvwr 2, to include variability in Fig1.
%follows the format of supp figure, with less panels, reduced size.
% showing ONLY the stride cycle variability across ppants.

    figure(1); clf;
set(gcf,'color', 'w', 'units', 'normalized', 'position', [.1 .1 .6, .69])
cd(procdatadir);

cd('GFX');

load('GFX_Data_inGaits', 'GFX_headY');
 % show all sub averages.
 fntsize=14
 subplot(2,3,2); cla
 allMin=[]
 for ippant = 1:size(GFX_headY,2);
     
     tmp =GFX_headY(ippant).doubgc_raw;
     tmp(tmp==0)=nan;
     
     % end plot at minimum point.
     [minIs, minAt] = min(tmp(80:end));
     minAt= minAt+79;
     plot(tmp(1:minAt), 'k');
     hold on
     allMin(ippant)=minAt;
 end
 ylim([-.06 .06]);
 title('Average per participant')
 ylabel({['detrended'];['head height (m)']});
% xlim([0 1.35])
prvX = get(gca,'xlim');
 set(gca,'xtick', 0:18:119,'xticklabel', [0:18:119]./90);
xlabel('Time (sec)')
set(gca,'fontsize', fntsize)
%
 subplot(2,3,3)
cfg.gap=.1;
% subplotBelow(gca,cfg)
histogram(allMin./90, 20); 
ylim([0 8])
text(.96,8*.9, ['(\itN\rm=' num2str(length(allMin)) ')'], 'fontsize',fntsize)
ylabel('Participant count')
xlabel('Stride duration (sec)')
set(gca,'fontsize', fntsize)
xlim([.95 1.35])
% xlim(prvX./90)
box on
 %%


end
 