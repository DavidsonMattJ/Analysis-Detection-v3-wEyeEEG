% j2C_RVW_append_doubGCData_shuffledVersion 
% building off the back of j2- > adding an extra column for the
% classification of an event in either a LRL or RLR step sequence. 
% necessary for showing double step proportions of DVs etc.

% append event percentages, relative to doubGait.


% now creating nPerm = 1000, where we randomly shuffle the stride-cycle
% onsets across the experiment.

% easiest way is to change the trial matching.

%laptop:
setmydirs_detectv3;
cd(procdatadir)
%% show ppant numbers:
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%
%%

% the analysis is too large, and will crash matlab. a solution is to create
% the surrogate data in one go.




for ippant= 1:nsubs
    cd(procdatadir)    %%load data from import job.
    load(pfols(ippant).name, ...
        'HeadPos', 'trialInfo', 'trial_summaryTable', 'subjID');
    savename = pfols(ippant).name;
    disp(['Preparing j2 cycle data... ' savename]);

    %% Gait extraction.
    % per event, find the relevant gait. then add the gPcnt relative
    % to LRL anr RLR sequence (if possible).

    allevents = size(trial_summaryTable,1);

    
    %%
    
    perm_strides_table=table();
    
    

% the key is to select the troughs from a different trial (at random).

newtrialOrder= randperm(100,100);
    for ievent= 1:allevents
        %guardclause to skip trial if  practice or stationary.
        if trial_summaryTable.isPrac(ievent) || trial_summaryTable.isStationary(ievent)
            %>>>>> add this new info to our table,
            %trgOnset:
%           
            continue
            perm_strides_table.(['trgO_gPcnt_LRL_p' num2str(iperm)])(ievent)= nan;
            perm_strides_table.(['trgO_gPcnt_RLR_p' num2str(iperm)])(ievent)= nan;
            
            
            %response
            perm_strides_table.(['respO_gPcnt_LRL_p' num2str(iperm)])(ievent)= nan;
            perm_strides_table.(['respO_gPcnt_RLR_p' num2str(iperm)])(ievent)= nan;


%             continue
        else
            % check if this trial is flagged for rejection (visual inspect
            % output of j1).
            
            skip=0;
            tmpitrial = trial_summaryTable.trial(ievent);
            itrial = newtrialOrder(tmpitrial);
            %store for debugging:
            perm_strides_table.observed_itrial(ievent) = tmpitrial;
            perm_strides_table.perm_itrial(ievent) = itrial;
            
            
            rejTrials_detectv3; %toggles skip based on bad trial ID
            if skip ==1
                continue;
            end

            %for each event, conver the target onset time, and reaction onset
            %time, to a % of a double gait cycle.

            % keep target onset times the same, but use different STRIDE
            % onsets.
            tTrial = trial_summaryTable.trial(ievent);

            % % ! ! !! now using the permuted troughs:
            
            trs = HeadPos(itrial).Y_gait_troughs;
            trs_sec = HeadPos(itrial).Y_gait_troughs_sec;
            pks = HeadPos(itrial).Y_gait_peaks;
            
            trialTime = trialInfo(itrial).times;
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
                    perm_strides_table.([savecols{ieventtype}  '_gPcnt_LRL_p' num2str(iperm)])(ievent)= nan;
                    perm_strides_table.([savecols{ieventtype}  '_gPcnt_RLR_p' num2str(iperm)])(ievent)= nan;
                   

                else % extract info we need:

                    % we have the gIndx, so look at prev and current step.
                    %
                    currentFt = trial_summaryTable.([savecols{ieventtype} '_gFoot'])(ievent); % note this was provided in j2B
                    if strcmp(currentFt, 'LR');
                        currDft = 'RLR';
                        prevDft = 'LRL';
                    else
                        currDft = 'LRL';
                        prevDft= 'RLR';
                    end

                    skipCurr=0;
                    skipPrev=0;
                    % extract the event as % double gait.[0 100]
                    
                    % might not have enough steps. 
                    if gindx>=2
                    prevSamps= trs(gindx-1):trs(gindx+1);
                    else 
                        prevSamps=nan;
                        skipPrev=1;
                    end
                    % same for current step
                    if gindx+2 <= length(trs)
                    currSamps= trs(gindx):trs(gindx+2);
                    else
                        currSamps=nan;
                        skipCurr=1;
                    end
%                     
%                     old way:
%                     resizeT= imresize(gaitTimes', [1,100]);
%                     gPcnt = dsearchn(resizeT', thisEvent);
                        useSamps = {currSamps, prevSamps};
                        doubgPcnts=nan(1,2);
                        for iwin=1:2
                            if skipCurr==1 && iwin==1
                                continue % pass over.
                            end
                            if skipPrev==1 && iwin==2
                                continue %pass aswell
                            end
                            gaitSamps= useSamps{iwin};
                            gaitTimes = trialInfo(tTrial).times(gaitSamps);
                        
                            % or take as proportion of total.
                            tDur = gaitTimes(end)- gaitTimes(1);
                            tE= thisEvent-gaitTimes(1);
                            gPcnt= round((tE/tDur)*100);
                            %
                            if gPcnt== 100
                                gPcnt=1;
                            end
                            doubgPcnts(iwin)= gPcnt;

                        end

                    %>>>>> add this new info to our table,
                    perm_strides_table.([savecols{ieventtype}  '_gPcnt_' currDft '_p' num2str(iperm)])(ievent)= doubgPcnts(1);
                    perm_strides_table.([savecols{ieventtype}  '_gPcnt_' prevDft '_p' num2str(iperm)])(ievent)= doubgPcnts(2);
                   


                end % if not misses
            end % targ and response onset

        end % if not practice

    end % each row in table (event)
    
% debug. Show trgOs
% ptrgOs = trial_summaryTable.trgO_gPcnt;
% ptrgOs = ptrgOs(~isnan(ptrgOs));
% hist(ptrgOs, 100);

%% include zscored version of RTs:
% allRTs = trial_summaryTable.clickRT;
%  userow = find(~isnan(allRTs));
% zRTs = zscore(allRTs(userow));
%    %place in table as new columnL
%    trial_summaryTable.z_clickRT = nan(size(trial_summaryTable,1),1);
%    trial_summaryTable.z_clickRT(userow) = zRTs;



   %% critical for later stages, remove the data for 'bad' trials.
   allbadtrials = badtrials; % pulled from rejTrials_detectv3.
   allts = trial_summaryTable.trial;
   remtrials = ismember(allts, allbadtrials);
   trial_summaryTable(remtrials,:)=[];
   
   end % for all permutations

%%
%%
   disp(['Finished appending DOUBLE gait percentage data for ... ' subjID]);
   save(savename, 'trial_summaryTable','-append');
end % participant

