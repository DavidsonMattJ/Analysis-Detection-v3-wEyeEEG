% j2_appendgaitcycleData

% Here we will load the summary data table, and add event markers as new columns, for when events
% happened relative to % gait cycle completion.
% For now, events of interest are target onset, and response onset. 

% jobs:
% -  append event percentages
% -  appends whether L/R ft (now in j2a).


%laptop:
setmydirs_detectv3;
cd(procdatadir)
%% show ppant numbers:
pfols = dir([pwd filesep 's0*']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%
%%


for ippant= 1:nsubs

    cd(procdatadir)
   load(pfols(ippant).name, ...
        'HeadPos', 'trialInfo', 'trial_summaryTable', 'subjID');
    savename = pfols(ippant).name;
    disp(['Preparing j2 cycle data... ' savename]);


    %% Gait extraction.
    % Per trial (and event), extract gait samples (trough to trough), normalize along x
    % axis, and store various metrics.


    allevents = size(trial_summaryTable,1);



    for ievent= 1:allevents
        %guardclause to skip trial if  practice or stationary.
        if trial_summaryTable.isPrac(ievent) || trial_summaryTable.isStationary(ievent)
            %>>>>> add this new info to our table,
            %trgOnset:
            trial_summaryTable.trgO_gCount(ievent) = nan;
            trial_summaryTable.trgO_gPcnt(ievent)= nan;
            trial_summaryTable.trgO_gDur(ievent)= nan;
             trial_summaryTable.trgO_gStart_sec(ievent)= nan;
              trial_summaryTable.trgO_gFin_sec(ievent)= nan;
               trial_summaryTable.trgO_gStart_samp(ievent)= nan;
                trial_summaryTable.trgO_gFin_samp(ievent)= nan;

            
            %response
            trial_summaryTable.respO_gCount(ievent) = nan;
            trial_summaryTable.respO_gPcnt(ievent)= nan;
            trial_summaryTable.respO_gDur(ievent)= nan;
               trial_summaryTable.respO_gStart_sec(ievent)= nan;
              trial_summaryTable.respO_gFin_sec(ievent)= nan;
               trial_summaryTable.respO_gStart_samp(ievent)= nan;
                trial_summaryTable.respO_gFin_samp(ievent)= nan;
      
              

            continue
        else
            % check if this trial is flagged for rejection (visual inspect
            % output of j1).
            skip=0;
            itrial = trial_summaryTable.trial(ievent);
            rejTrials_detectv3; %toggles skip based on bad trial ID
            if skip ==1
                continue;
            end

            %for each event, conver the target onset time, and reaction onset
            %time, to a % of a gait cycle.

            tTrial = trial_summaryTable.trial(ievent);

            trs = HeadPos(tTrial).Y_gait_troughs;
            trs_sec = HeadPos(tTrial).Y_gait_troughs_sec;
            pks = HeadPos(tTrial).Y_gait_peaks;
            trialTime = trialInfo(tTrial).times;
            eventcolumns = {'targOnset', 'targRT'};
            savecols = {'trgO', 'respO'};

            for ieventtype=1:2
                %%Target onsets first:
                % find the appropriate gait. Final trough that the event is later
                % than:
                thisEvent = trial_summaryTable.([ eventcolumns{ieventtype} ])(ievent);
                
               gindx = find(thisEvent>trs_sec,1,'last');
                %don't fill for misses, or before/after first/last step in a trial.
                if thisEvent ==0 || isempty(gindx) || gindx== length(trs)  %
                    trial_summaryTable.([savecols{ieventtype}  '_gCount'])(ievent) =nan;
                    trial_summaryTable.([savecols{ieventtype}  '_gPcnt'])(ievent)= nan;
                    trial_summaryTable.([savecols{ieventtype}  '_gDur'])(ievent)= nan;
                    trial_summaryTable.([savecols{ieventtype}  '_gStart_sec'])(ievent) = nan;
                    trial_summaryTable.([savecols{ieventtype}  '_gFin_sec'])(ievent) = nan;
                    trial_summaryTable.([savecols{ieventtype}  '_gStart_samp'])(ievent) = nan;
                    trial_summaryTable.([savecols{ieventtype}  '_gFin_samp']) (ievent)= nan;

                else % extract info we need:

                    if isnan(thisEvent)
                        pause
                    end
                    % extract the event as % gait.[0 100]
                    gaitsamps =trs(gindx):trs(gindx+1); % avoid counting edge bins twice.
                    gaitTimes = trialInfo(tTrial).times(gaitsamps);
                    resizeT= imresize(gaitTimes', [1,100]);
                    gPcnt = dsearchn(resizeT', thisEvent);

                    % or take as proportion of total. 
                    tDur = gaitTimes(end)- gaitTimes(1);
                    tE= thisEvent-gaitTimes(1);
                    gPcnt= round((tE/tDur)*100);

                    %>>>>> add this new info to our table,
                    trial_summaryTable.([savecols{ieventtype}  '_gCount'])(ievent) = gindx; %trgO_gcount, or respO_gcount.
                    trial_summaryTable.([savecols{ieventtype}  '_gPcnt'])(ievent)= gPcnt;
                    trial_summaryTable.([savecols{ieventtype}  '_gDur'])(ievent)= round(gaitTimes(end)-gaitTimes(1),3);
                    trial_summaryTable.([savecols{ieventtype}  '_gStart_sec'])(ievent) = [gaitTimes(1)];
                    trial_summaryTable.([savecols{ieventtype}  '_gFin_sec'])(ievent) = [gaitTimes(end)];
                      trial_summaryTable.([savecols{ieventtype}  '_gStart_samp'])(ievent) = [gaitsamps(1)];
                    trial_summaryTable.([savecols{ieventtype}  '_gFin_samp'])(ievent) = [gaitsamps(end)];


                end % if not misses
            end % targ and response onset

        end % if not practice

    end % each row in table (event)
 
 

 
% debug. Show trgOs
% ptrgOs = trial_summaryTable.trgO_gPcnt;
% ptrgOs = ptrgOs(~isnan(ptrgOs));
% hist(ptrgOs, 100);

%% include zscored version of RTs:
allRTs = trial_summaryTable.clickRT;
 userow = find(~isnan(allRTs));
zRTs = zscore(allRTs(userow));
   %place in table as new columnL
   trial_summaryTable.z_clickRT = nan(size(trial_summaryTable,1),1);
   trial_summaryTable.z_clickRT(userow) = zRTs;



   %% critical for later stages, remove the data for 'bad' trials.
   allbadtrials = badtrials; % pulled from rejTrials_detectv3.
   allts = trial_summaryTable.trial;
   remtrials = ismember(allts, allbadtrials);
   trial_summaryTable(remtrials,:)=[];
%%
%%
   disp(['Finished appending gait percentage data for ... ' subjID]);
   save(savename, 'trial_summaryTable','-append');
end % participant

