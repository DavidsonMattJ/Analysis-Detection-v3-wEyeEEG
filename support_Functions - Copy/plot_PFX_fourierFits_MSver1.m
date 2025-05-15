
function     plot_PFX_fourierFits_MSver1(cfg, testData, PFX_FourierNull);
% plots a 2 x 3 figure (probably supplementary), showing the individual
% fits overlayed, and a histogram for peak per participant.

%%
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


%set up figure:
figure(1); clf;
set(gcf, 'color', 'w', 'units','normalized', 'position', [.1 .1 .6 .6])
hDataAll=[];
for id=1:3

%which field of datastructure to plot?
if strcmp(cfg.DV{id}, 'RT')
    usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_rts'];
elseif strcmp(cfg.DV{id}, 'Accuracy')
    usefield = [gaitfield{nGaits_toPlot}  binfield{usebin+1} '_Acc'];

elseif strcmp(cfg.DV{id}, 'Counts')
    usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_counts'];
end

ppantData=[];

for ippant = 1:length(PFX_FourierNull)

ppantData(ippant,:) = PFX_FourierNull(ippant).([cfg.type{id} 'Ons_' usefield '_fitsRsq_Obs']);
end % ippant


subplot(3,3,pc);
% plot overlyed version:
    plot(cfg.Hzspace, ppantData', 'color',[.8 .8 .8])
hold on;
plot(cfg.Hzspace, mean(ppantData,1), barCols{id}, 'LineWidth',2)
    %h
    
    ylabel('R^2');
    xlabel('Frequency (cycles-per-gait)');
set(gca,'fontsize',fntsize)
% per participant, find the max in Hzspace. plot a histogrma.
[mx, idx]= max(ppantData,[],2);
hdata = cfg.Hzspace(idx);
% subplot(2,3,4:6); hold on


subplot(3,3,pc+3)
imagesc(cfg.Hzspace, 1:size(ppantData,1), ppantData)
xlabel('Frequency (cycles-per-gait)');
ylabel('participant id');
c=colorbar;
ylabel(c, 'R^2');
clim([0 .4])
% hg=histogram(hdata, 0:.5:10);
% hg.FaceColor= barCols{id};
% hg.FaceAlpha= .2;
%  ylabel('Participant count');
%  xlabel('R^2 peak (Frequency)');
pc=pc+1;
% hDataAll(id,:)= histcounts(hdata, 0:.5:10);
% ylim([0 11])
set(gca,'fontsize',fntsize)

end % id

% subplot(2,3,4:6); hold on; cla
% hg=bar(.5:.5:10, hDataAll, 'stacked');


end% function