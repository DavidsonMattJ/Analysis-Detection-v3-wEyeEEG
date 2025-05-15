
function plot_PFX_prevalence_MSver2(cfg, testdata, PFX_FourierNull)
% % % within this function, we will display the population prevalence.
% at present, only accuracy data is available (permutations at ppant level
% are slow).

nGaits_toPlot=2;
pc=1; % plot counter
barCols= {'b', 'r', 'm'};
% different colours to distinguish from the DVs 

% blues:
barCols= {[127,0,255]./255, [56,109,249]./255, [19, 200,231]./255,[91,249,200]./255, [237,200,11]./255} ; % blues
%sunset
barCols = {[113,29,176]./255, [194,18,146]./255, [239, 64,64]./255,[255,167,50]./255};


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

figure(1); clf
set(gcf, 'color', 'w', 'units','normalized', 'position', [.05 .05 .65 .95])
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

if rankBy2Hz
% only rank once.
% if id==1

at2 = dsearchn(cfg.Hzspace', 2);
[i,nOrder]=sort(ppantData(:, at2));
% end

ppantData= ppantData(nOrder,:);
ppantNHST= ppantNHST(nOrder,:)
end

% let's colour by whether 2 or 4 Hz was stronger?
compHz_2cps= dsearchn(cfg.Hzspace', [1.5,2.5]');
compHz_4cps= dsearchn(cfg.Hzspace', [3.5,4.5]');
compHz_5cps= dsearchn(cfg.Hzspace', [5]');

pCols= {};
%%
% clf
[count2, count4, countOther, countNone]=deal(0);

for ippant = 1:size(ppantData,1)


 
 if any(ppantNHST(ippant,compHz_2cps(1):compHz_2cps(2)))==1 % sig at 2 hz.
         subplot(5,3,1 + (id-1));
         colp= barCols{1};
         count2 = count2+1;
 elseif any(ppantNHST(ippant,compHz_4cps(1):compHz_4cps(2)))==1 
         subplot(5,3,4 + (id-1));
     colp= barCols{2};
     count4= count4+1;
 elseif any(ppantNHST(ippant,compHz_5cps(1):end))==1
     subplot(5,3,7 + (id-1));
         colp= barCols{3};
         countOther=countOther+1;
 elseif sum(ppantNHST(ippant, :))==0 % no sig
  
       subplot(5,3,7 + (id-1));
         colp= [.2 .2 .2];
         countNone=countNone+1;
     
 end
 hold on 

p=plot(cfg.Hzspace, ppantData(ippant,:),'color',colp, 'LineWidth',2);
p.Color = [p.Color, .7];
ylim([0 .7]);

end
% add counts:
subplot(5,3,1 + (id-1));
text(4, .75, ['\bf\itn\rm\bf=' num2str(count2) ' \rm(NHST at ~2 cps)' ],'color', barCols{1});
set(gca,'xtick', round(cfg.Hzspace(1:5:length(cfg.Hzspace)),1), 'fontsize',12);
ylabel('R^2');
subplot(5,3,4 + (id-1));
text(4, .75,  ['\bf\itn\rm\bf=' num2str(count4) ' \rm(NHST at ~4 cps)' ], 'color', barCols{2});
set(gca,'xtick', round(cfg.Hzspace(1:5:length(cfg.Hzspace)),1), 'fontsize',12);
ylabel('R^2');
subplot(5,3,7 + (id-1));
set(gca,'xtick', round(cfg.Hzspace(1:5:length(cfg.Hzspace)),1), 'fontsize',12);
xlabel('Frequency (cps)');
ylabel('R^2');
text(4, .75,  ['\bf\itn\rm\bf=' num2str(countOther) ' \rm(NHST > 5 cps)' ], 'color', barCols{3});
text(4, .65, ['\bf\itn\rm\bf=' num2str(countNone) ' \rm(no oscillations)']);
set(gca,'fontsize', 12)
%%

%%
% axis ij

% subplot(412);
%attach to base

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

%% summarise some results:
% check NHST at 2, and 4, 6 Hz.
checkme= dsearchn(cfg.Hzspace', [2,4,6]');
someF= sum(ppantNHST(:, checkme),2); 

allSig = length(find(someF));
nMulti = length(find(someF>1));
% any sig oscillation:
anyF  = sum(ppantNHST,2); 
N_anySig = length(find(anyF));
%%
% c
%%

%compute bayes prev on total (any freq sig).
% now plot the prevalence results 
% for 2, 4, other, and any.
nPrevEsts = [count2, count4, N_anySig];

useCols= {barCols{1}, barCols{2}, 'k'};% {'b', 'r', 'k'};
nPrevcycs= {'2', '4', 'any'};
useSts= {'-', '-', '-'};
% pl
legh=[];
% plot the dummy version:
if id==1
subplot(5,3, [(10 +(id-1)), (13 +(id-1))]) ; hold on

posterior = bayesprev_posterior(x,count2,n,alpha);
% add MAP as a point
xmap = bayesprev_map(count2,n,alpha);
pmap = bayesprev_posterior(xmap,count2,n,alpha);
bayesprev_hpdi(0.95,count2,n,alpha);
% plot(x, posterior,'Color',useCols{iPrev}, 'linew',2, 'linest',useSts{iPrev} );
pMAP=plot(.45, 9,'.','MarkerSize',30,'Color',[.7 .7 .7]);
hold on
oil = 2; % linewidths. (outer)
iil = 4; % (iNNER)
pHDPI= plot([.4 .5],[9 9],'Color',[.7 .7 .7],'LineWidth',oil)
hold on;
text(.4, 9, '      MAP [96% HDPI]', 'EdgeColor','k', 'fontsize',10)
end 
for iPrev = 1:length(nPrevEsts);

    k = nPrevEsts(iPrev);
posterior = bayesprev_posterior(x,k,n,alpha);
% add MAP as a point
xmap = bayesprev_map(k,n,alpha);
pmap = bayesprev_posterior(xmap,k,n,alpha);
h95 = bayesprev_hpdi(0.95,k,n,alpha);

% now plot!
% figure(2); 
%
subplot(5,3, [(10 +(id-1)), (13 +(id-1))]) ; hold on
% subplot(1,3,id); hold on
legh(iPrev)=plot(x, posterior,'Color',useCols{iPrev}, 'linew',2, 'linest',useSts{iPrev} );
plot(xmap, pmap,'.','MarkerSize',30,'Color',useCols{iPrev});

plot([h95(1) h95(2)],[pmap pmap],'Color',useCols{iPrev},'LineWidth',oil)

% add 50% HPDI
% h50 = bayesprev_hpdi(0.5,k,n,alpha);
% plot([h50(1) h50(2)],[pmap pmap],'Color',useCols{iPrev},'LineWidth',iil)

xlabel('Population prevalence proportion')
ylabel('Posterior density');
ylim([0 10])
textmsg= ['     ' sprintf('%.2f', xmap) ' [' sprintf('%.2f', h95(1)) ', ' sprintf('%.2f', h95(2)) ']' ];
textmsg=['   Sig. at '  nPrevcycs{iPrev} ' cps']
if iPrev==3
    textmsg = ['   \itany\rm oscillations'];
end

% txt=text(.05, 10-iPrev, textmsg, 'EdgeColor', 'k');
% plot(.075, 10-iPrev,'.','MarkerSize',30,'Color',useCols{iPrev});
set(gca,'fontsize',12)

if iPrev==3
text(xmap, pmap+1, [sprintf('%.2f', xmap) ' [' sprintf('%.2f', h95(1)) ', ' sprintf('%.2f', h95(2)) ']' ], ...
    'HorizontalAlignment', 'center');

end
end
%%
lgd=legend([legh],{'2 cps','4 cps','\itAny'},'Location','NorthWest');
% lgd.EntryContainer.IconDisplayStyle = 'off';  % Hide the default referenced objects
lgd.ItemTokenSize(1) = 8;
%%
% %
% disp([' subjects with more than one sig freq: ' num2str(nMulti)])
% disp([' subjects with ANY sig freq: ' num2str(allSig)])
% set(gca,'fontsize', cfg.fntsize);
%%


GFX_NHST(id,:,:) =ppantNHST;

end % id
%% determine overlap.
% GFX_NHST (DV, ppants, freqs).
GFX_DVcount=[];
bandsAre={compHz_2cps, compHz_4cps};
for ippant = 1:size(GFX_NHST,2) 
    
    % see whether we have the same freq across DVs.
    for iDV=1:3
        for iband= 1:2
         
    GFX_DVcount(ippant,iDV,iband) = any(squeeze(GFX_NHST(iDV,ippant,bandsAre{iband}(1):bandsAre{iband}(2))));
        end
    end

end
% how many had sig in more than one band, per DV:
DVcount = squeeze(sum(GFX_DVcount,3)); % look for indexs >1
nMultiperAcc = length(find(DVcount(:,1)>1));
nMultiperRT = length(find(DVcount(:,2)>1));
nMultiperResp = length(find(DVcount(:,3)>1));
%
disp('Count of ppants with sig at multiple bands for :')
disp(['Accuracy data : n = ' num2str(nMultiperAcc)]);
disp(['RT data : n = ' num2str(nMultiperRT)]);
disp(['Resp data : n = ' num2str(nMultiperResp)]);

%coutns of participants with oscillations at different bands across DVs:
Hzcount = squeeze(sum(GFX_DVcount,2)); % now we have 2 columns, showing the nDVs sig at 2, or 4 cps.
% look for the 1,1 case (indicates no repetitions).
% find by replacing 0 with nan, and look for sum=1;
Hzcount(Hzcount==0)=nan;
NMixedHz = length(find(sum(Hzcount,2)==2));

disp('Count of ppants with different sig CPS for different DVs :')
disp(['n= ' num2str(NMixedHz)]);



end