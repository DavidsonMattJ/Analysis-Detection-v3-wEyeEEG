
function plot_PFX_prevalence_MSver1(cfg, testdata, PFX_FourierNull)
% % % within this function, we will display the population prevalence.
% at present, only accuracy data is available (permutations at ppant level
% are slow).

nGaits_toPlot=2;
pc=1; % plot counter
barCols= {'b', 'r', 'm'};

% both this and the next use the same figure function:
iLR=3;
gaitfield = {'gc', 'doubgc'};
gaitprint = {'gait', 'stride'};
binfield = {'','_binned'};
usebin=1;
legp=[]; % for legend
fntsize=cfg.fntsize;

rankBy2Hz= 1; % rank the erp image by R2 strength.
%%
%set up figure:

figure(1); clf;
set(gcf, 'color', 'w', 'units','normalized', 'position', [.05 .05 .6 .95])
tsAre= {'Accuracy', 'Reaction time', 'Response onset'};
GFX_NHST=[]; % use to find overlap
%%
for id=1:3

%which field of datastructure to plot?
if strcmp(cfg.DV{id}, 'RT')
    usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_rts'];
elseif strcmp(cfg.DV{id}, 'Accuracy')
    usefield = [gaitfield{nGaits_toPlot}  binfield{usebin+1} '_Acc'];

elseif strcmp(cfg.DV{id}, 'Counts')
    usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_counts'];
end
%%
ppantData=[];
ppantBounds=[];
ppantNHST=[];
for ippant = 1:length(PFX_FourierNull)

ppantData(ippant,:) = PFX_FourierNull(ippant).([cfg.type{id} 'Ons_' usefield '_fitsRsq_Obs']);

ppantBounds(ippant,:) =PFX_FourierNull(ippant).([cfg.type{id} 'Ons_' usefield '_fitsRsq_ShuffCV'])(3,:); % 3rd pos is 95%

% binary NHST:
ppantNHST(ippant,:) = ppantData(ippant,:) > ppantBounds(ippant,:);
end % ippant
%%
%% plot all bounds.
figure(1); 
nppants= size(ppantData,1);
subplot(3,3,id)
if rankBy2Hz
% only rank once.
% if id==1

at2 = dsearchn(cfg.Hzspace', 2);
[i,nOrder]=sort(ppantData(:, at2));
% end

ppantData= ppantData(nOrder,:);
end
imagesc(cfg.Hzspace, 1:size(ppantData,1), ppantData);
%
caxis([0 .35]);

c=colorbar;
%
ylabel(c, 'R^2', 'units', 'normalized', 'Position', [c.Position(1)+.45, c.Position(2)+.05], ...
  'Rotation',0, 'VerticalAlignment','top', 'HorizontalAlignment','center'); shg

set(gca,'xtick', [], 'fontsize', cfg.fntsize)

%%
% axis ij
% xlabel('Frequency (cpg)')
ylabel('Participant');
set(gca,'ytick', 0:4:32, 'yticklabel', num2str(nOrder), 'ydir', 'normal')
title(tsAre{id}, 'HorizontalAlignment', 'right', 'fontsize',cfg.fntsize);
shg
%

% subplot(412);
%attach to base
gax= gca;
offsetY=.01;
%%
% new position directly beneath prev plot:
oldpos = gax.Position; % x, y, width, height
newpos = [oldpos(1), (oldpos(2)-oldpos(4)/2-offsetY), oldpos(3),  (oldpos(4)/2 - offsetY)];
subplot('position', newpos);
%%
% subplot(3,3,id+3)

plot(cfg.Hzspace, mean(ppantNHST,1), 'o-k');
ylabel({['NHST'];['Prevalence']}, 'fontsize', cfg.fntsize);

xlabel('Frequency (cps)', 'fontsize', cfg.fntsize);
set(gca,'xtick', round(cfg.Hzspace(1:4:length(cfg.Hzspace)),1))
shg
%%
% perform preavalence estimate.
% we have the first level, participant NHST result (0 or 1 at each
% frequency).

% Bayesian prevalence inference is performed with three numbers: 
% k, the number of significant participants (e.g. sum of binary indicator
% variable)
% n, the number of participants in the sample
% alpha, the false positive rate
NSub = size(ppantNHST,1);


bayes_postResults=[];
sigref=[];
for ifreq = 1:length(cfg.Hzspace)
   indsig = ppantNHST(:,ifreq);
k = sum(indsig);
n = NSub;
alpha = 0.05; % default value see 'help ttest'

% plot posterior distribution of population prevalence

co = get(gca,'ColorOrder'); 
ci=1;
hold on

x = linspace(0,1,100);
posterior = bayesprev_posterior(x,k,n,alpha);
% plot(x, posterior,'Color','k');

% add MAP as a point
xmap = bayesprev_map(k,n,alpha);
pmap = bayesprev_posterior(xmap,k,n,alpha);
% plot(xmap, pmap,'.','MarkerSize',20,'Color',co(ci,:));

% add lower bound as a vertical line
bound = bayesprev_bound(0.95,k,n,alpha);
boundY= bayesprev_posterior(bound,k,n,alpha);
% line([bound bound], [0 bayesprev_posterior(bound,k,n,alpha)],'Color',co(ci,:),'LineStyle',':')

% add 95% HPDI
% oil = 2; % linewidths. (outer)
% iil = 4; % (iNNER)
h95 = bayesprev_hpdi(0.95,k,n,alpha);
% plot([h95(1) h95(2)],[pmap pmap],'Color',co(ci,:),'LineWidth',oil)
% add 50% HPDI
h50 = bayesprev_hpdi(0.5,k,n,alpha);
% plot([h(1) h(2)],[pmap pmap],'Color',co(ci,:),'LineWidth',iil)
% subplot()

% xlabel('Population prevalence proportion')
% ylabel('Posterior density')

%>>Store.
bayes_postResults(ifreq).posterior = posterior;
bayes_postResults(ifreq).xmap = xmap; % prevalence MAP
bayes_postResults(ifreq).pmap = pmap; % prevalence posterior
bayes_postResults(ifreq).bound= bound;
bayes_postResults(ifreq).boundY= boundY;

bayes_postResults(ifreq).h95= h95; % prevalence posterior
bayes_postResults(ifreq).h50= h50; % prevalence posterior
bayes_postResults(ifreq).Sig = h95(1)>0;
sigref(ifreq)= h95(1)>0;

end %
%% add those significant to the map: 
if id==1 % accuracy cols
rmap= cbrewer('seq', 'Reds', 10);
bmap= cbrewer('seq', 'Blues', 10);
cmap = [rmap(5:9,:);bmap(5:10,:)];
elseif id==2
rmap= cbrewer('seq', 'Reds', 10);
bmap= cbrewer('seq', 'Blues', 10);
pmap= cbrewer('seq', 'Purples', 10);
gmap = cbrewer('seq', 'Greens',10);

cmap = [rmap(5:9,:);bmap(5:9,:); pmap(5:7,:); gmap(7,:)];
elseif id==3
rmap= cbrewer('seq', 'Reds', 10);
bmap= cbrewer('seq', 'Blues', 10);
pmap= cbrewer('seq', 'Purples', 10);
gmap = cbrewer('seq', 'Greens',10);

cmap = [rmap(5:9,:);bmap(5:8,:); pmap(5:7,:)];



end
x = linspace(0,1,100);


psig = find(sigref);

% backtrack (add to previous plot).
for isig= 1:length(psig)
p=plot(cfg.Hzspace(psig(isig)), mean(ppantNHST(:, psig(isig)),1), 'o', 'color', cmap(isig,:), 'MarkerFaceColor',cmap(isig,:));
% plot(cfg.Hzspace(psig(isig)), mean(ppantNHST(:, psig(isig)),1), 'x', 'color', cmap(isig,:));

end

% next figures (prevalence stats:)

postleg=[];

%create location for next plot
subplot(6,3,[10,13] +(id-1));
% gax= gca;
% offsetY=.005;
%%
% new position directly beneath prev plot:
% oldpos = gax.Position; % x, y, width, height
% newpos = [oldpos(1), (oldpos(2)-oldpos(4)/2-offsetY), oldpos(3),  (oldpos(4) + offsetY+ oldpos(4)/2)];
% subplot('position', newpos);


% first plot a dummy version for the legend:
   xmap = bayes_postResults(psig(1)).xmap;
   pl=  plot(xmap, 1,'.','MarkerSize',40,'Color',[.5 .5 .5]);
   h95 = bayes_postResults(psig(1)).h95;
hold on;
   pr= plot([h95(1) h95(2)],[1 1],'Color',[.5 .5 .5],'LineWidth',2);


   % now overplot:
for isig= 1:length(psig)
    indx= psig(isig);
    thisHz = cfg.Hzspace(indx);

      xmap = bayes_postResults(indx).xmap;
hold on
    plot(xmap, isig,'.','MarkerSize',40,'Color',cmap(isig,:));

h95 = bayes_postResults(indx).h95;

    plot([h95(1) h95(2)],[isig isig],'Color',[cmap(isig,:), .5],'LineWidth',2)
% add 50% HPDI
h50 = bayes_postResults(indx).h50;
plot([h50(1) h50(2)],[isig isig],'Color',[ cmap(isig,:), .5],'LineWidth',4)

% text(.75, isig, [sprintf('%.2f',xmap) ',[' sprintf('%.2f', h95(1)) ', ' sprintf('%.2f', h95(2)) ']']);
% xlabel('Population prevalence proportion')
xlim([0 1]);
ylim([cfg.Hzspace(1), cfg.Hzspace(end)]);
ylabel('Frequency (cps)','fontsize', cfg.fntsize)
set(gca,'ydir', 'reverse', 'fontsize', cfg.fntsize)
end
%%
set(gca, 'Ytick', 1:(length(psig)),'YTickLabel', string(round(cfg.Hzspace(psig),1)))
ylim([0 length(psig)+1])

if id==1
 legend([pl, pr], {'MAP', '95% HDPI'});
end
%  set(gca, 'Ytick', 1:(length(psig)),'YTickLabel', string(cfg.Hzspace(psig)))


shg
%%
cd(cfg.figdir);
cd('Population prevalence');

print('-dpng', [cfg.DV{id} ' prevalence']);

GFX_NHST(id,:,:) =ppantNHST;

%% summarise some results:
% check NHST at 2, and 4, 6 Hz.
checkme= dsearchn(cfg.Hzspace', [2,4,6]');
anyF= sum(ppantNHST(:, checkme),2); 

allSig = length(find(anyF));
nMulti = length(find(anyF>1));
%%
% cla
% add to plot?
txt1= ['N_a_n_y = ' num2str(allSig) '/' num2str(n) ];
txt2= ['N_m_u_l_t_i = ' num2str(nMulti) '/' num2str(allSig)]; 
txt=text(.65, length(psig)-1, txt1, 'Interpreter','tex');
txt=text(.65, length(psig), txt2, 'Interpreter','tex');
shg
%%

%compute bayes prev on total (any freq sig).
k = sum(allSig);


posterior = bayesprev_posterior(x,k,n,alpha);
% add MAP as a point
xmap = bayesprev_map(k,n,alpha);
pmap = bayesprev_posterior(xmap,k,n,alpha);
h95 = bayesprev_hpdi(0.95,k,n,alpha);

% now plot!
% figure(2); 
%%
subplot(6,3, 16 +(id-1)); hold on
% subplot(1,3,id); hold on
plot(x, posterior,'Color','k');
plot(xmap, pmap,'.','MarkerSize',40,'Color',barCols{id});

oil = 2; % linewidths. (outer)
iil = 4; % (iNNER)
plot([h95(1) h95(2)],[pmap pmap],'Color',barCols{id},'LineWidth',oil)

% add 50% HPDI
h50 = bayesprev_hpdi(0.5,k,n,alpha);
plot([h50(1) h50(2)],[pmap pmap],'Color',barCols{id},'LineWidth',iil)

xlabel('Population prevalence proportion')
ylabel('Posterior density');
ylim([0 7])
%
textmsg= [sprintf('%.2f', xmap) ' [' sprintf('%.2f', h95(1)) ', ' sprintf('%.2f', h95(2)) ']' ];
% add details.
% textbox = annotation('textbox', [0.2, 0.7, 0.2, 0.2], 'String', textmsg, 'EdgeColor', 'k', 'LineWidth', 1.5);

txt=text(.05, 4, textmsg, 'EdgeColor', 'k');
disp([' subjects with more than one sig freq: ' num2str(nMulti)])
disp([' subjects with ANY sig freq: ' num2str(allSig)])
set(gca,'fontsize', cfg.fntsize);
%%
end % id
%% determine overlap.
% GFX_NHST (DV, ppants, freqs).
GFX_DVcount=[];
for ippant = 1:size(GFX_NHST,2) 
    
    % see whether we have the same freq across DVs.
    
    DVatFreq = sum(squeeze(GFX_NHST(:,ippant,checkme)),1);
    % this is now the number of DVs that were sig at each freq.
    GFX_DVcount(ippant,:) = DVatFreq;
    
end
% just eyeballing (in command window). The number of instances when there
% is only a single freq sig on a DV, and a different (singular) freq,
%is present for another DV, is n=6. 
%e.g. [ 0 0 1] = one DV sig 6 cps 
%e.g. [2 0 0 ] = 2 DVs sig 2 cps.
%e.g. [0 1 1] = 1 DV sig at 4 cps, 1 DV sig at 6 cps **

% N only one DV oscillation:
Nonly1 = length(find(sum(GFX_DVcount,2) ==1));
Nonly1 = length(find(sum(GFX_DVcount,2) ==1));


end