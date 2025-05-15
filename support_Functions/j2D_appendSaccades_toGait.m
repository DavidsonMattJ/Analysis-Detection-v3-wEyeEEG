% j2D_appendSaccades_toGait
% job2D_appendSaccades_toGait


% load the eye movement data per participant, and identify saccade onset
% (and duration).
% create a table and assign when saccades happened relative to the
% co-occurring gait.
% note that saccades are ID'd in previous script (job1B_timestampSaccades)

%gait start/fin was performed in j2_epochHead_byssteps.

%output is a table with saccades per trial, with gait percent .
%while here we could also look at saccade direction toward/away from
%targets ons creen.

setmydirs_detectv3;

cd(procdatadir);

visualiseResults=0;
subjFols= dir([pwd filesep '*summary_data.mat']);
%%
for ippant=1:length(subjFols)
  cd(procdatadir)
 
    savename= subjFols(ippant).name;
    
    %new table per ppant.
    saccadeSummaryTable=table();
    saccadeTS_raw=nan(1,3,100); % this will be a largematrix, raw time series sacc start and end (X,Y dims).
    saccadeTS_relative=nan(1,3,100); % this will be a largematrix, time series from 0,0 at sacc start .

    iSaccadeCounter= 1; % this is the row in our table.


        % load the head pos and eye pos data (prev imported).
        load([savename], 'EyeDir', 'HeadPos', 'EyeOrigin', 'trial_summaryTable');
         
        %trialsummarytable contains the target onset times (and types)
        trialSummaryTable= trial_summaryTable; % back convert
        ftmp = find(savename =='_');
        subjID = savename(1:ftmp(1)-1);

        disp(['Preparing ' subjID]);
    
       
        for itrial= 1:length(HeadPos)
        
            % skip if identified as a poor trial:
            skip=0;
            rejTrials_detectv3;
            if skip ==1
                continue;
            end


            % let's be conservative, and apply a minimum saccade duration
            % of 2 samples (22 ms).
            saccX = find(EyeDir(itrial).saccade_dursX>=2);
            saccY = find(EyeDir(itrial).saccade_dursY>=2);

            SaccadeEvents =[saccX, saccY];
          
            if isempty(SaccadeEvents)
                continue
            end
            if isempty(HeadPos(itrial).Y_gait_troughs)
                walkSpeed=0;
                continue
            else
                walkSpeed=1;

            end
            
            % if we have saccades, add each to our table:
            % add to table with relevant gait information:
           
            
            %for each event, determine whether the gait cycle is Left-Right
            %ft, or Right-Left ft. This is done by observing the sway (left
            % to right), in head position.
            
           
            % Head position data:
            tmpPos=  squeeze(HeadPos(itrial).Y);
            tmpSway = squeeze(HeadPos(itrial).Z);
            tmpwalkDir = squeeze(HeadPos(itrial).X);

            % quick classification:
            if mod(itrial,2)~=0 % odd numbers (walking toward computer)
                % Then descending on X, (i.e. first trajectory), more positive z values
                % are left side of the body.
                walkDir= 'IN';
                ylab='Right - > Left';
            elseif mod(itrial,2) ==0
                % ascending (the return trajectory), pos z values are RHS
                walkDir= 'OUT';

                ylab='Left - > Right';
            end
             % is the z value increasing or decreasing, relative to feet placement?

                
                %%
               

            trs = HeadPos(itrial).Y_gait_troughs;
%             trs_sec = HeadPos(itrial).Y_gait_troughs_sec;
            pks = HeadPos(itrial).Y_gait_peaks;

            % lets add saccadeOnset, offset(landing)
            % and saccade RT.
            flds1= {'Y_clean', 'Z_clean'};
            flds2= {'saccade_startsY', 'saccade_startsX'};
            flds3 = {'saccade_dursY', 'saccade_dursX'};
            %
            saccIndex= {saccY, saccX};
            savecols= {'Y', 'X'};
            for idim=1:2 % Y then X

                sOnsets= EyeDir(itrial).([flds2{idim}]);
                sDurs= EyeDir(itrial).([flds3{idim}]);
                sEnds = sOnsets+sDurs;
                %convert to seconds:                
%                 sStarts = timevec(sOnsets);
%                 sDurs= timevec(sDurs);
%                 sEnds = sStarts+sDurs;


                saccEvents= saccIndex{idim};
                
                % we can also get some target information:
                relRows = find(trialSummaryTable.trial==itrial);

                targOnsets = trialSummaryTable.targOnset(relRows);  

                for ievent = 1:length(saccEvents)

            thisOnset =  sOnsets(saccEvents(ievent));
            thisDur =  sDurs(saccEvents(ievent));
            thisOffset =  sEnds(saccEvents(ievent));
        
            % the above are all in samps.

            % find nearest gait
            gindx_Onset = find(thisOnset>trs,1,'last');
            gindx_Offset = find(thisOffset>trs,1,'last');
        

            % CRUCIAL! Note that eye movement data is recorded in world
            % coords. You can se this, as the z data is either +1 or -1
            % based on the direction of travel. (relative coords would
            % always be +1 facing forward). 
            %Therefore we need to adjust the sign of X locations based on
            %direction of travel.
            % if z= 1, lets assume positive x is RHS of fixation
%             if median(EyeDir(itrial).X_clean)>0 % positive, no change.
% 
%             else % participant has turned 180, a negative target location is always on the LHS
%                 % but X dir of eye pos needs to be updated.
%                 EyeDir(itrial).Z_clean(thisOnset)=EyeDir(itrial).Z_clean(thisOnset)*-1; 
%             end

            % also save the starting and ending locations (in cartesian space).
            saccStart_loc = [EyeDir(itrial).X_clean(thisOnset), EyeDir(itrial).Y_clean(thisOnset),EyeDir(itrial).Z_clean(thisOnset)];
            saccEnd_loc = [EyeDir(itrial).X_clean(thisOffset), EyeDir(itrial).Y_clean(thisOffset),EyeDir(itrial).Z_clean(thisOffset)];
            
            
            % centre saccade at 0-0.
           
            saccEnd_locO = saccEnd_loc - (saccStart_loc);
            
            % extract TS while we are here.
            saccadeTimeseries=[];
            saccadeTimeseries(1,:) = EyeDir(itrial).X_clean(thisOnset:thisOffset);
            saccadeTimeseries(2,:) = EyeDir(itrial).Y_clean(thisOnset:thisOffset);
            saccadeTimeseries(3,:) = EyeDir(itrial).Z_clean(thisOnset:thisOffset);
             adjTS=[];
            adjTS(1,:) = saccadeTimeseries(1,:) - saccadeTimeseries(1,1);
            adjTS(2,:) = saccadeTimeseries(2,:) - saccadeTimeseries(2,1);
            adjTS(3,:) = saccadeTimeseries(3,:) - saccadeTimeseries(3,1);


           

            %%

            %don't fill for misses, or before/after first/last step in a trial.
                if isempty(gindx_Onset) || gindx_Onset== length(trs)  %
                    if iSaccadeCounter>1 % suppress auto warnings.
                         saccadeSummaryTable(iSaccadeCounter,:)= saccadeSummaryTable(iSaccadeCounter-1,:);
                    end
                    
                    saccadeSummaryTable.(['walkSpeed'])(iSaccadeCounter) =walkSpeed;
                       saccadeSummaryTable.(['trialID'])(iSaccadeCounter) =itrial;
                       saccadeSummaryTable.(['saccadeOnsetFrame'])(iSaccadeCounter) =nan;
                       saccadeSummaryTable.('saccOns_gFoot')(iSaccadeCounter) = {nan};
                       saccadeSummaryTable.('saccOffs_gFoot')(iSaccadeCounter) = {nan};

                        saccadeSummaryTable.(['gindx_Onset'])(iSaccadeCounter) =nan;  
                        saccadeSummaryTable.(['gindx_Offset'])(iSaccadeCounter) =nan;  
                      % which dim was it detected on?
                    saccadeSummaryTable.(['saccDim']) (iSaccadeCounter)= {savecols{idim}};
                    saccadeSummaryTable.(['saccOnset_gPcnt']) (iSaccadeCounter)= nan;
                     saccadeSummaryTable.(['saccOffset_gPcnt']) (iSaccadeCounter)= nan;
                      saccadeSummaryTable.(['saccDur_samps']) (iSaccadeCounter)= nan;
                    saccadeSummaryTable.(['saccStart_LocX']) (iSaccadeCounter)= nan;
                    saccadeSummaryTable.(['saccStart_LocY']) (iSaccadeCounter)= nan;
                    saccadeSummaryTable.(['saccStart_LocZ']) (iSaccadeCounter)= nan;
                    saccadeSummaryTable.(['saccEnd_LocX']) (iSaccadeCounter)= nan;
                    saccadeSummaryTable.(['saccEnd_LocY']) (iSaccadeCounter)= nan;
                    saccadeSummaryTable.(['saccEnd_LocZ']) (iSaccadeCounter)= nan;

%                        saccadeSummaryTable.(['saccStart_Loc']) (iSaccadeCounter)= nan;
                      saccadeSummaryTable.(['saccEnd_LocX_relO']) (iSaccadeCounter)= nan;
                      saccadeSummaryTable.(['saccEnd_LocY_relO']) (iSaccadeCounter)= nan;
                      saccadeSummaryTable.(['saccEnd_LocZ_relO']) (iSaccadeCounter)= nan;
%                       saccadeSummaryTable.(['saccEnd_LocX_relTrg']) (iSaccadeCounter)= nan;
%                       saccadeSummaryTable.(['saccEnd_LocY_relTrg']) (iSaccadeCounter)= nan;
%                       saccadeSummaryTable.(['saccEnd_LocZ_relTrg']) (iSaccadeCounter)= nan;
%                       
%                       saccadeSummaryTable.(['currentTarget_idx']) (iSaccadeCounter)= nan;
%                       saccadeSummaryTable.(['saccTowardTarg']) (iSaccadeCounter)= nan;
%                       saccadeSummaryTable.(['saccadeRT'])(iSaccadeCounter)=nan;
%                       saccadeSummaryTable.(['saccadeRT2'])(iSaccadeCounter)=nan;
%                        saccadeSummaryTable.(['saccadeRT3'])(iSaccadeCounter)=nan;
%                   
                      % will remove later
                      saccadeTS_raw(iSaccadeCounter,:,1:length(saccadeTimeseries)) = nan;
                      saccadeTS_relative(iSaccadeCounter,:,1:length(saccadeTimeseries)) = nan;
                else % extract info we need:

                   
                    % extract the event as % gait.[0 100]
                   %
                   useIndx= [gindx_Onset, gindx_Offset];
                   useSamp = [thisOnset, thisOffset];
                   for iOnsOff=1:2
                   
                       % sacconset:
                       gindx = useIndx(iOnsOff);
                       if gindx~= length(trs)
                           gaitsamps =trs(gindx):trs(gindx+1); % avoid counting edge bins twice.
                           tDur = gaitsamps(end)- gaitsamps(1);
                           tE= useSamp(iOnsOff)-gaitsamps(1);
                           gPcnt_= round((tE/tDur)*100);
                       else
                           gPcnt_= nan;
                       end


                    % here we can also append which foot they were
                    % swinging:
                    %
                    gStart=gaitsamps(1);
                    gEnd= gaitsamps(end);
                midlineS = mean(tmpSway([gStart,gEnd]));
                meanS = mean(tmpSway(gStart:gEnd));

                gaitSway = tmpSway(gStart:gEnd);

                % here we can also append which foot they were
                % swinging:
                %

            if strcmp(walkDir, 'IN') && meanS>midlineS
                % if we were approaching computer, and swinging positive
                ft='LR'; % foot placement was right-left.
            elseif strcmp(walkDir, 'IN') && meanS<midlineS
                %swinging opp direction.
                ft= 'RL';
            elseif strcmp(walkDir, 'OUT') && meanS>midlineS
                %if back to computer, positive swing is LR
                ft= 'RL';
            elseif strcmp(walkDir, 'OUT') && meanS<midlineS
                ft= 'LR';
            end
            


                % rename for correct save:
                if iOnsOff==1
                    gPcnt_Onset= gPcnt_;
                    ft_Onset = ft;
                else
                    gPcnt_Offset= gPcnt_;
                    ft_Offset = ft;
                end
            

            end % iOns Off
                

                %% sanity check:
%                 clf; subplot(211);
% 
%                 plot(1:length(tmpSway), tmpSway);
%                 subplot(212);
%                 plot(1:length(tmpSway), tmpSway, 'color', [.8 .8 .8]);
%                 % to ease interp, reorient to allocentric:
%                 if strcmp(walkDir, 'OUT')
%                     set(gca,'Ydir','reverse')
%                      ylabel('Right -Left')
%                 else
%                     ylabel(ylab)
%                 end
% % 
%                 %overlay
%                 hold on;
%                 plot(gStart:gEnd, gaitSway, 'k', 'linew', 2)
%                 hold on;
%                 plot([gStart, gEnd], tmpSway([gStart, gEnd]), 'ro')
%                 %add text to be sure
%                 title(['Walking ' walkDir ', gait: ' ft])
%                



                    % note that for each saccade, we have a start, end
                    % location, we could compute if it is toward / away
                    % from the target? maybe another script.

                    %>>>>> add this new info to our table,
                    if iSaccadeCounter>1 % suppress auto warnings.
                        saccadeSummaryTable(iSaccadeCounter,:)= saccadeSummaryTable(iSaccadeCounter-1,:);
                    end
                           saccadeSummaryTable.(['walkSpeed'])(iSaccadeCounter) =walkSpeed;
                       saccadeSummaryTable.(['trialID'])(iSaccadeCounter) =itrial; 
                       saccadeSummaryTable.(['saccadeOnsetFrame'])(iSaccadeCounter) =thisOnset; 

                       saccadeSummaryTable.('saccOns_gFoot')(iSaccadeCounter) = {ft_Onset};
                       saccadeSummaryTable.('saccOffs_gFoot')(iSaccadeCounter) = {ft_Offset};

                        saccadeSummaryTable.(['gindx_Onset'])(iSaccadeCounter) =gindx_Onset;  
                        saccadeSummaryTable.(['gindx_Offset'])(iSaccadeCounter) =gindx_Offset;  
                      
                    saccadeSummaryTable.(['saccDim']) (iSaccadeCounter)= {savecols{idim}};
                    saccadeSummaryTable.(['saccOnset_gPcnt']) (iSaccadeCounter)= gPcnt_Onset;
                     saccadeSummaryTable.(['saccOffset_gPcnt']) (iSaccadeCounter)= gPcnt_Offset;
                      saccadeSummaryTable.(['saccDur_samps']) (iSaccadeCounter)= thisDur;
 
                      % cartesian coordinates (raw)
                      saccadeSummaryTable.(['saccStart_LocX']) (iSaccadeCounter)= saccStart_loc(1);
                    saccadeSummaryTable.(['saccStart_LocY']) (iSaccadeCounter)= saccStart_loc(2);                    
                    saccadeSummaryTable.(['saccStart_LocZ']) (iSaccadeCounter)= saccStart_loc(3);
                    
                    saccadeSummaryTable.(['saccEnd_LocX']) (iSaccadeCounter)= saccEnd_loc(1);
                    saccadeSummaryTable.(['saccEnd_LocY']) (iSaccadeCounter)= saccEnd_loc(2);
                    saccadeSummaryTable.(['saccEnd_LocZ']) (iSaccadeCounter)= saccEnd_loc(3);

                    % diff between sacc end and sacc start:
                    saccadeSummaryTable.(['saccEnd_LocX_relO']) (iSaccadeCounter)= saccEnd_locO(1);
                      saccadeSummaryTable.(['saccEnd_LocY_relO']) (iSaccadeCounter)= saccEnd_locO(2);
                      saccadeSummaryTable.(['saccEnd_LocZ_relO']) (iSaccadeCounter)= saccEnd_locO(3); % this bugs show trying:
%                       saccadeSummaryTable.(['saccEnd_LocZ_relO']) (iSaccadeCounter)= saccStart_locO(3); % this bugs show trying:

                      %diff between sacc end and target loc:
%                       saccadeSummaryTable.(['saccEnd_LocX_relTrg']) (iSaccadeCounter)= saccEnd_locTrg(1);
%                       saccadeSummaryTable.(['saccEnd_LocY_relTrg']) (iSaccadeCounter)= saccEnd_locTrg(2);
%                       saccadeSummaryTable.(['saccEnd_LocZ_relTrg']) (iSaccadeCounter)= saccEnd_locTrg(3);
%                       
% 
%                       saccadeSummaryTable.('currentTarget_idx') (iSaccadeCounter)= thistargLoc;
%                       saccadeSummaryTable.('saccTowardTarg') (iSaccadeCounter)= saccTowardTarg; % 0 or 1
%                       saccadeSummaryTable.('saccadeRT')(iSaccadeCounter)=saccRT; % onset since most recent target,
%                        saccadeSummaryTable.('saccadeRT2')(iSaccadeCounter)=saccRT2; % onset since n-1 target,
%                         saccadeSummaryTable.('saccadeRT2')(iSaccadeCounter)=saccRT3; % onset since n-1 target,
%                     % will remove later
                      saccadeTS_raw(iSaccadeCounter,:,1:length(saccadeTimeseries)) = saccadeTimeseries; % all saccades made
                      saccadeTS_relative(iSaccadeCounter,:,1:length(saccadeTimeseries)) = adjTS;

                end % if appropriate gait
                % 
                iSaccadeCounter=iSaccadeCounter+1;
                % prefill next line to avoid warning output in command
                % window.
               
                end % each event in saccade dimension
            end % X Y dim
            
% %% visualize saccadeS?
% if visualiseResults==1
% clf;
% flds1= {'Y_clean', 'Z_clean'};
% flds2= {'saccade_startsY', 'saccade_startsX'};
% flds3 = {'saccade_dursY', 'saccade_dursX'};
% cols= {'b', 'r'};
% for idim= 1:2
%     subplot(2,1,idim);
%     plot(timevec, EyeDir(itrial).([flds1{idim}]), cols{idim}, 'LineWidth',2); hold on;
% 
%     sOnsets= EyeDir(itrial).([flds2{idim}]);
%     sDurs= EyeDir(itrial).([flds3{idim}]);
% 
%     sStarts = timevec(sOnsets);
%     sDurs= timevec(sDurs);
%     sEnds = sStarts+sDurs;
%     yyaxis right
%     for iS = 1:length(sOnsets)
%         % highlight on the plot the start and duration.
% 
%         plot([sStarts(iS), sStarts(iS)], [-.5 .5], [cols{idim} '-o'])
% 
%         % add a patch?
% 
%         xP = [sStarts(iS) sStarts(iS) sEnds(iS) sEnds(iS)];
%         yP=[-.5 .5 .5 -.5];
% 
%         ph= patch(xP, yP, [1 1 1], 'FaceColor', cols{idim});
%         ph.FaceAlpha=.2;
%             end
% end
% end
%%

      

        end % all trials

 %% critical for later stages, remove the data for 'bad' trials.
   allbadtrials = badtrials; % pulled from rejTrials_detectv3.
   allts = saccadeSummaryTable.trialID;
   remtrials = ismember(allts, allbadtrials);
   saccadeSummaryTable(remtrials,:)=[];
   saccadeTS_relative(remtrials,:,:)=[];
   saccadeTS_raw(remtrials,:,:)=[];

 
%% tidy up by removing nan trials:
nanrows = find(isnan(saccadeSummaryTable.saccadeOnsetFrame));
saccadeSummaryTable(nanrows,:)=[];
saccadeTS_relative(nanrows,:,:)=[];
saccadeTS_raw(nanrows,:,:)=[];

% sanity checks
% clf;
if visualiseResults==1
figure(1); clf
uset = find(~isnan(saccadeSummaryTable.saccOnset_gPcnt));
nbins=50;
subplot(311)
histogram(saccadeSummaryTable.saccOnset_gPcnt(uset),nbins); title('sacc onset')
subplot(312)
histogram(saccadeSummaryTable.saccOffset_gPcnt(uset),nbins); title('sacc offset')

% for the third, plot all RTs by gIndx.
RTpergIndx= nan(1,100);
for ipcnt=1:100
    subtrsA = find(saccadeSummaryTable.saccOnset_gPcnt==ipcnt);
    subtrsB = find(saccadeSummaryTable.saccTowardTarg==1);
      subtrs= intersect(subtrsA, subtrsB);
    meanRT = nanmean(saccadeSummaryTable.saccadeRT(subtrs));
    RTpergIndx(ipcnt)= meanRT;
end
subplot(313)
bar(1:100, RTpergIndx); title('sacc RT')


%%
figure(2);
clf
for itrial=1:length(uset)
    pltrial= uset(itrial);
    subplot(121);
    trialD= squeeze(saccadeTS_raw(pltrial,:,:));
    plot(trialD(3,:), trialD(2,:), 'color', [.5 .5 .5 .5]); hold on
title('raw saccade direction')
    subplot(122);
    trialD= squeeze(saccadeTS_relative(pltrial,:,:));
    plot(trialD(3,:), trialD(2,:), 'color', [.5 .5 .5 .5]); hold on
    title('relative saccade direction')

end
%
  stimPos=[-.2 .2; .2 .2 ; .2 -.2; -.2 -.2; ...
            -.07 .07; .07 .07; .07 -.07; -.07 -.07 ];
%%
figure(3); clf
dirCols= {'r', 'g'}; % away toward
for iloc=1:8
loctrials = find(saccadeSummaryTable.currentTarget_idx==iloc);

subtrials = intersect(loctrials,uset);
% further subdivide.?
for iT=1:2
    subss = find(saccadeSummaryTable.saccTowardTarg==(iT-1));
    subtrials_byDir = intersect(subtrials,subss);

 subplot(2,4,iloc);
 plot(stimPos(iloc,1), stimPos(iloc,2), 'bd', 'LineWidth',4);
    for itrial= 1:length(subtrials_byDir)
       
        saccEndX = saccadeSummaryTable.saccEnd_LocX_relO(subtrials_byDir(itrial));
        saccEndY = saccadeSummaryTable.saccEnd_LocY_relO(subtrials_byDir(itrial));

%         saccEndX = saccadeSummaryTable.saccEnd_LocX(subtrials_byDir(itrial));
%         saccEndY = saccadeSummaryTable.saccEnd_LocY(subtrials_byDir(itrial));


%         saccEndX = saccadeSummaryTable.saccEnd_LocX_relTrg(subtrials_byDir(itrial));
%         saccEndY = saccadeSummaryTable.saccEnd_LocY_relTrg(subtrials_byDir(itrial));
        
        plot(saccEndX, saccEndY,'o', 'color', dirCols{iT}); hold on
    end
end
hold on; plot([0 0], ylim, 'k-')
hold on; plot(xlim, [0 0], 'k-')
    xlim([-1 1]); ylim([-1 1])
%     hold on;
%     plot(squeeze(mean(saccadeTS_raw(subtrials,1,:),1)),...
%         squeeze(mean(saccadeTS_raw(subtrials,2,:),1)), 'r-')

 subplot(2,4,iloc);
 plot(stimPos(iloc,1), stimPos(iloc,2), 'bd', 'LineWidth',4);
end
end
%%
save([savename],  'saccadeSummaryTable','-append');

disp(['Finished j2D for ' subjID]);
end % participant
