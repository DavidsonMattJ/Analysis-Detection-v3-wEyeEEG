
function plot_PFX_prevalence(cfg, testdata, PFX_FourierNull)
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
usebin=1
legp=[]; % for legend
fntsize=12;


rankBy2Hz= 1; % rank the erp image by R2 strength.
%%

for id=3%:3
    %set up figure:
figure(1); clf;
set(gcf, 'color', 'w', 'units','normalized', 'position', [.1 .1 .6 .6])


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
subplot(2,3,1)
if rankBy2Hz

at2 = dsearchn(cfg.Hzspace', 2);
[i,n]=sort(ppantData(:, at2));

ppantData= ppantData(n,:);
end
imagesc(cfg.Hzspace, 1:size(ppantData,1), ppantData);
%
caxis([0 .35]);

c=colorbar;
ylabel(c, 'R^2', 'Position', [c.Position(1)+c.Position(3), c.Position(2)+.15], ...
  'Rotation',0, 'VerticalAlignment','top', 'HorizontalAlignment','center'); shg
set(gca,'xtick', [])
%%
% axis ij
% xlabel('Frequency (cpg)')
ylabel('Participant');
set(gca,'ytick', 0:4:32, 'ydir', 'normal')
shg
%

% subplot(412);
%attach to base
gax= gca;
offsetY=.005
%%
% new position directly beneath prev plot:
oldpos = gax.Position; % x, y, width, height
newpos = [oldpos(1), (oldpos(2)-oldpos(4)/2-offsetY), oldpos(3),  (oldpos(4)/2 - offsetY)];
subplot('position', newpos);
%%
% subplot(3,3,id+3)

plot(cfg.Hzspace, mean(ppantNHST,1), 'o-k');
ylabel({['NHST'];['Prevalence']});

xlabel('Frequency (cpg)');
set(gca,'xtick', cfg.Hzspace(1:4:length(cfg.Hzspace)))
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


for isig= 1:length(psig)
    indx= psig(isig);
    thisHz = cfg.Hzspace(indx);

    % add posterior distribution:
subplot(2,3,2)
   plot(x, bayes_postResults(indx).posterior, 'color',cmap(isig,:) ); hold on;
      hold on; 
      % add lowerbound
      bound = bayes_postResults(indx).bound;
     line([bound bound], [0 bayes_postResults(indx).boundY],'Color',cmap(isig,:) ,'LineStyle',':')

      % add MAP point
      xmap = bayes_postResults(indx).xmap;
      pmap =  bayes_postResults(indx).pmap;

   postleg(isig)=  plot(xmap, pmap,'.','MarkerSize',20,'Color',cmap(isig,:));

% xlabel('Population prevalence proportion')
ylabel('Posterior density')

gax= gca;
offsetY=.005;
%%
% new position directly beneath prev plot:
oldpos = gax.Position; % x, y, width, height
newpos = [oldpos(1), (oldpos(2)-oldpos(4)/2-offsetY), oldpos(3),  (oldpos(4)/2 - offsetY)];
subplot('position', newpos);

hold on
    plot(xmap, isig,'.','MarkerSize',20,'Color',cmap(isig,:));

h95 = bayes_postResults(indx).h95;

    plot([h95(1) h95(2)],[isig isig],'Color',cmap(isig,:),'LineWidth',2)
% add 50% HPDI
h50 = bayes_postResults(indx).h50;
plot([h50(1) h50(2)],[isig isig],'Color',cmap(isig,:),'LineWidth',4)

xlabel('Population prevalence proportion')
xlim([0 1]);
ylabel('Frequency (cpg)')
set(gca,'ydir', 'reverse')
end

set(gca, 'Ytick', 1:length(psig),'YTickLabel', string(cfg.Hzspace(psig)))
subplot(2,3,2)
 legend(postleg, string(cfg.Hzspace(psig)), 'fontsize',8)

shg



end % id






end