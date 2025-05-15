% j4B_crunchPFX_concatGFX_relativetoSaccades


% Here we load all data from previous analyses, preparing to save before
% plotting in next round of jobs.

%%%%%% QUEST DETECT version 3 (wEEG and eye) %%%%%%



setmydirs_detectv3;
% show ppant numbers:
cd(procdatadir)
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%   

job.concatGFX=1;
job.plotGFX=1;

%% concat data:
%preallocate storage:

pidx1=ceil(linspace(1,100,21)); % length n-1
pidx2=ceil(linspace(1,100,41));% 
useindexes = {pidx1, pidx2};
%

%first load the presaved GFX_SaccadeData
GFX_TargRelativeToSaccadePos=[];
subjIDs=cell([1, length(pfols)]);

if job.concatGFX==1
for ippant =1:nsubs
    cd(procdatadir)    
    load(pfols(ippant).name, 'subjID', 'saccadeSummaryTable','trial_summaryTable','HeadPos');
    
    subjIDs{ippant} = subjID;
    
    disp(['concatenating subject saccade relative data ' subjID]);
    
    
  %% NOW  we will process the saccade onset data (where saccades occurred during the gait). 
    % also calculate binned versions
    % 
    ppantSaccData= saccadeSummaryTable;
    ppantGaitData= trial_summaryTable;
    
    % pseudoCode. 
    %for each target onset in our table, add a column denoting position
    %reltive to saccade onset.

    for ievent = 1:size(trial_summaryTable,1 )

        %find nearest saccade in same trial, 
        thisOnset = trial_summaryTable.targOnset(ievent);
        thistrial=trial_summaryTable.trial(ievent);
        thisBlock = trial_summaryTable.block(ievent); 

        % find all the saccades that occurred within this trial, and the
        % nearest for labelling:
        saccTrialIndx =find(saccadeSummaryTable.trialID==thistrial);
    
        if isempty(saccTrialIndx)
            %no saccades on this trial for alignment.
            trial_summaryTable.trgO_preSacc(ievent)=nan;
            trial_summaryTable.trgO_postSacc(ievent)=nan;

        else
        
            % find the nearest locs.
            saccOnsets_sec = saccadeSummaryTable.saccadeOnsetFrame(saccTrialIndx)./90;
            % convert to sec, and find nearest.
            % last saccade before this onset (for post sacc window).
            postSacc= find(saccOnsets_sec<thisOnset,1,'last');
            preSacc= find(saccOnsets_sec>thisOnset,1,'first'); % pre sacc
            
            if isempty(postSacc)
                trial_summaryTable.trgO_postSacc(ievent)=nan;
            else
                trgO_postsacc= thisOnset- saccOnsets_sec(postSacc);

                % to make things easy, round to 2dp;
                trial_summaryTable.trgO_postSacc(ievent)=round(trgO_postsacc,2);
            end

            if isempty(preSacc)
                trial_summaryTable.trgO_preSacc(ievent)=nan;
            else
                trgO_presacc= thisOnset- saccOnsets_sec(preSacc);
                trial_summaryTable.trgO_preSacc(ievent)=round(trgO_presacc,2);
            end


        end



    end
            
    disp(['fin append for ippant ' num2str(ippant)])

%% now we have to work out the average per bin. 
uniqueBins = unique([trial_summaryTable.trgO_preSacc; trial_summaryTable.trgO_postSacc]);
uniqueBins(isnan(uniqueBins))=[];


posIdx = [trial_summaryTable.trgO_preSacc; trial_summaryTable.trgO_postSacc];
allRTs = [trial_summaryTable.clickRT;trial_summaryTable.clickRT];
allCor = [trial_summaryTable.targCor;trial_summaryTable.targCor];

% PFX_trgrelativetoSacc=[];
binCounts = [];
binRTs=[];
binAcc=[];

for ibin = 1:length(uniqueBins)

    usebin = uniqueBins(ibin);

    allEvents = find(posIdx==usebin);
    
    binCounts(ibin)= length(allEvents);
    binRTs(ibin)= nanmean(allRTs(allEvents));
    binAcc(ibin)= nanmean(allCor(allEvents));

end % bin
% 
%%
% clf; figure(1)
% xvec = uniqueBins;
% subplot(311);
% bar(xvec, binCounts);
% xlim([-4 4])
% % xvec = uniqueBins;
% subplot(312);
% bar(xvec, binRTs);
% xlim([-4 4])
% subplot(313);
% bar(xvec, binAcc);
% xlim([-4 4])
% shg
%%
% truncate to store

binRange = [-4:.01:4];
% due to missing data within some bin indices, we need to store
% individually:
[PFX_counts, PFX_RTs, PFX_Acc]=deal(nan(1,length(binRange)));

% for storage just use +-4 sec.
binLims = dsearchn(uniqueBins,[-4 4]');

for ibin = binLims(1):binLims(2)
    
    storeBin = uniqueBins(ibin);
    storeIdx = dsearchn(binRange', storeBin);
    PFX_counts(storeIdx) = binCounts(ibin);
    PFX_RTs(storeIdx) = binRTs(ibin);
    PFX_Acc(storeIdx)  = binAcc(ibin);

end

%% check:
% figure(2); clf
% xvec = binRange;
% subplot(311);
% bar(xvec, PFX_counts);
% 
% % xvec = uniqueBins;
% subplot(312);
% bar(xvec, PFX_RTs);
% 
% subplot(313);
% bar(xvec, PFX_Acc);
% shg
%%




GFX_TargRelativeToSaccadePos.saccbinRange = binRange;
GFX_TargRelativeToSaccadePos.saccbinCounts(ippant,:) = PFX_counts;
GFX_TargRelativeToSaccadePos.saccbinRTs(ippant,:) =PFX_RTs;
GFX_TargRelativeToSaccadePos.saccbinAcc(ippant,:) =PFX_Acc;

end % ppant


%%
%% save GFX
cd([procdatadir filesep 'GFX']);
disp('Saving GFX');
save('GFX_Data_inGaits', ...
   'GFX_TargRelativeToSaccadePos',...
   '-append');%, '-append');

end 


if job.plotGFX==1
    cd(procdatadir); cd('GFX'); load('GFX_Data_inGaits.mat', 'GFX_TargRelativeToSaccadePos');

%% clf
figure(1);clf; set(gcf,'color','w','units','normalized','Position',[.1 .1 .8 .8])


% ploteach field
fldsare= {'saccbinCounts', 'saccbinRTs', 'saccbinAcc'};
ysare= {'Target Counts', 'Target reaction time [s]', 'Hit Rate (%)'};

colsAre= {'k', 'r', 'b' };
for ifield=1:3

    plotData = GFX_TargRelativeToSaccadePos.([fldsare{ifield}]);

    plotM =squeeze(nanmean(plotData,1));
    plotErr = CousineauSEM(plotData);

%     plotM= smooth(plotM)
    subplot(1,3,ifield);
    bh=bar(GFX_TargRelativeToSaccadePos.saccbinRange, plotM);
    bh.FaceColor = colsAre{ifield};
    hold on;

    errorbar(GFX_TargRelativeToSaccadePos.saccbinRange, plotM, plotErr,'LineStyle','none')
xlim([-.5 .5]);
xlabel('Time from saccade onset (sec)');
ylabel(ysare{ifield});
set(gca,'fontsize',15);
hold on;
plot([0 0 ], ylim, 'k:')
end
shg



end

