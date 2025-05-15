% AA_TargetDistribution_inTrial\

%%%%%% QUEST & EEG version

% Detection experiment (contrast)
%%  Import from csv. FramebyFrame, then summary data.

%laptop:
cd('C:\GitHub\Analysis-Detection-v3-wEyeEEG\support_Functions');
%
sourcedir = pwd;
cd ../Figures
figdir= pwd;
cd ../data_Processed;
procdatadir = pwd;

%%
% quick count and average for group stats.

GFX_diffTOns=[];
GFX_TOnsHist=[];
GFX_TargCounts=[];
GFX_TargCounts_sum=[];
GFX_TargCounts_stat=[];
GFX_TargCounts_walk=[];
for ippant=1:length(pfols)
    % cd(datadir)

    pfols = dir([pwd filesep '*summary_data.mat']);

    %% load subject data as table.
    filename = pfols(ippant).name;
    %extract name&date from filename:
    ftmp = find(filename =='_');
    subjID = filename(1:ftmp(end)-1);


    savename = [subjID '_summary_data'];

    load(filename, 'trial_summaryTable', 'trialInfo');

    %% plot walk:

    % for each walk type (7 or 8). plot distribution of trial target onsets.
    alltrials = unique(trial_summaryTable.trial);
    all7 =[];
    all8=[];
    allT= [];
    diffTs=[];
    allCounts =[];
    allCounts_Walk=[];
    allCounts_Stationary=[];
    minIs=0;
    maxIs=1;
    for itrial = 1:length(alltrials)

        user= find(trial_summaryTable.trial == itrial);
        tOns = trial_summaryTable.targOnset(user);
        %
        %     if length(tOns)==7
        %         all7 = [all7, tOns'];
        %     elseif length(tOns)>7
        %         all8 = [all8, tOns'];
        %     end

        allT= [allT, tOns'];
        diffTs= [diffTs, diff(tOns)']; % average difference between targets.
        allCounts= [allCounts, length(user)];

        if ~isempty(user)
            if trial_summaryTable.isStationary(user(1)) ==1
                allCounts_Stationary = [allCounts_Stationary, length(user)];
            else
                allCounts_Walk = [allCounts_Walk, length(user)];
            end

        end
    end

    if min(allCounts)<minIs
        minIs= min(allCounts);
    end
    if max(allCounts)> maxIs
        maxIs = max(allCounts);
    end
    %
    % clf;
    % subplot(311);
    % histogram(all7,[0:.05:9.5])
    % subplot(312)
    % histogram(all8,[0:.05:9.5] );
    % subplot(313)
    % histogram(allT,[0:.05:9.5] );

    allTcounts=histcounts(allT, [0:.05:9.5]);

    GFX_diffTOns(ippant)=mean(diffTs);
    GFX_TOnsHist(ippant,:)= allTcounts;
    GFX_TargCounts(ippant)= mean(allCounts);

    % if trial_summaryTable.isStationary(user(1))
    GFX_TargCounts_sum(ippant)= sum(allCounts);
    GFX_TargCounts_walk(ippant)= sum(allCounts_Walk);
    GFX_TargCounts_stat(ippant)= sum(allCounts_Stationary);
end
%
disp(mean(GFX_TargCounts_sum))
%%
% plot summary
clf;
bar(mean(GFX_TOnsHist,1));





