% j3_epochGait_timeseries

% using the trough information in the summary table (or HeadPos), we will
% now epoch the raw, resampled, and normalized versions of the gait
% time-series (for later plots).



% now running AFTER we have cleaned eye movement data, to epoch at the same
% time.
% ensure J1A_EyeMovementpertrial.m has been completed.

setmydirs_detectv3;
cd(procdatadir)
%% show ppant numbers:
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%
 

resampSize=100; 


for ippant =1:nsubs
    cd(procdatadir)    %%load data from import job.
    load(pfols(ippant).name, ...
        'HeadPos', 'trialInfo', 'trial_summaryTable', 'subjID', ...
        'EyePos', 'EyeDir'); % % new additions.
    savename = pfols(ippant).name;
    disp(['Preparing j3 gait time series data... ' savename]);

    %% Gait extraction.
    % Per trial (and event), extract gait samples (trough to trough), normalize along x
    % axis, and store various metrics.


    allevents = size(trial_summaryTable,1);

    %create a large empty matrix to staw raw gait, and resampled sizes
    %also:

    gait_ts_raw= zeros(1,500); % in samples
    [gait_ts_resamp, gait_ts_eyedirY_resamp,...
        gait_ts_eyedirZ_resamp,gait_ts_eyeposY_resamp,...
        gait_ts_eyeposZ_resamp] = deal(zeros(1,100)); % resampled to 100%


    % lets also grab the double GC:
    doubgait_ts_raw= zeros(1,500);
    doubgait_ts_resamp= zeros(1,200);

    [trialallocation, gaitDuration,gaitIdx, gaitSamps, gaitStart, gaitFin]= deal(nan);
    gaitFeet= {'L'};
    gait_ts_gData= table(trialallocation,gaitIdx, gaitFeet,  gaitStart,gaitSamps,gaitDuration);
    
%find all gait onsets (in samples)
gOnsets = find(~isnan(trial_summaryTable.trgO_gStart_samp));

% % This works for single Gaits only:
    for igait = 1:length(gOnsets)

        rowIndx = gOnsets(igait);
        itrial = trial_summaryTable.trial(rowIndx);
        
        %check if this trial should be skipped:
        skip=0;
        rejTrials_detectv3; %toggles skip based on bad trial ID (itrial)
        if skip ==1
            continue;
        end
        % continue:
        gaitStFin = [trial_summaryTable.trgO_gStart_samp(rowIndx),trial_summaryTable.trgO_gFin_samp(rowIndx)];
        gaitDur = trial_summaryTable.trgO_gFin_sec(rowIndx)- trial_summaryTable.trgO_gStart_sec(rowIndx);
        gFt = trial_summaryTable.trgO_gFoot(rowIndx);

        % % % now using a detrended version (smoother plots).
        % for each gait we have, extract the raw head time series,
       

% detrend the gait to avoid slant effects?
headY = detrend(HeadPos(itrial).Y);
% headY = (HeadPos(itrial).Y);
rawHead = headY(gaitStFin(1):gaitStFin(2));
Eye_Y = detrend(EyePos(itrial).Y_clean);
rawEye_Y= Eye_Y(gaitStFin(1):gaitStFin(2));

% while here, also extract eye origin and gaze direction.
%  rawHead = HeadPos(itrial).Y(gaitStFin(1):gaitStFin(2));        
rawGaze_Y = EyeDir(itrial).Y_clean(gaitStFin(1):gaitStFin(2));
rawGaze_Z = EyeDir(itrial).Z_clean(gaitStFin(1):gaitStFin(2));
% rawEye_Y = EyePos(itrial).Y_clean(gaitStFin(1):gaitStFin(2));
rawEye_Z = EyePos(itrial).Z_clean(gaitStFin(1):gaitStFin(2));

        %also resample along x dimension:
        resampHead = imresize(rawHead', [1,resampSize]);
        resampGY= imresize(rawGaze_Y', [1,resampSize]);
        resampGZ= imresize(rawGaze_Z', [1,resampSize]);
        resampOY= imresize(rawEye_Y', [1,resampSize]);
        resampOZ= imresize(rawEye_Z', [1,resampSize]);

        
        
        %gaitcount(within trial)
        gIdx = find(HeadPos(itrial).Y_gait_troughs==(gaitStFin(1)));
        
        %>>>>> store:
        gait_ts_raw(igait,1:length(rawHead)) = rawHead;
        gait_ts_resamp(igait,:) = resampHead;
        gait_ts_eyedirY_resamp(igait,:)= resampGY;
        gait_ts_eyedirZ_resamp(igait,:)= resampGZ;
        gait_ts_eyeposY_resamp(igait,:)= resampOY;
        gait_ts_eyeposZ_resamp(igait,:)= resampOZ;

        gait_ts_gData.trialallocation(igait) = itrial;
        gait_ts_gData.gaitDuration(igait)= gaitDur;
        gait_ts_gData.gaitFeet(igait) = gFt;
        gait_ts_gData.gaitIdx(igait) = gIdx;
        gait_ts_gData.gaitStart(igait) = gaitStFin(1);
        gait_ts_gData.gaitSamps(igait) = length(rawHead);
        %prefill next row to suppress output warnings in command window:
        gait_ts_gData(igait+1,:) = gait_ts_gData(igait,:) ;
    end
    %remove last row (autofilled)
    gait_ts_gData(igait+1,:)=[];
    %convert zeros to nan to tidy plot:
    gait_ts_raw(gait_ts_raw==0)=nan;
    gait_ts_resamp(gait_ts_resamp==0)=nan;

%% % now also epoch the doubgc Head data, by iterating over all gaits and epoching n+1.



  for igait = 1:length(gOnsets)

        rowIndx = gOnsets(igait);
        itrial = trial_summaryTable.trial(rowIndx);
        
        %check if this trial should be skipped:
        skip=0;
        rejTrials_detectv3; %toggles skip based on bad trial ID (itrial)
        if skip ==1
            continue;
        end
        % continue by looking for gait +1.

        gaitIdx = trial_summaryTable.trgO_gCount(rowIndx);
        allGs = HeadPos(itrial).Y_gait_troughs;
        %select either the current or prev started step for this doubGC.
        try
            gaitStFin= allGs(gaitIdx):allGs(gaitIdx+2);
        catch
            gaitStFin= allGs(gaitIdx-1):allGs(gaitIdx+1);
        end
            
%      

% detrend the gait to avoid slant effects?
headY = detrend(HeadPos(itrial).Y);
% headY = (HeadPos(itrial).Y);
rawHead = headY(gaitStFin);

        %also resample along x dimension:
        resampHead = imresize(rawHead', [1,200]); % double the length.
 %>>>>> store:
        doubgait_ts_raw(igait,1:length(rawHead)) = rawHead;
        doubgait_ts_resamp(igait,:) = resampHead;





  end
    %convert zeros to nan to tidy plot:
    doubgait_ts_raw(gait_ts_raw==0)=nan;
    doubgait_ts_resamp(gait_ts_resamp==0)=nan;












disp(['saving gait and double gait time series data for ' subjID])
cd(procdatadir);
save(savename, 'gait_ts_raw', 'gait_ts_resamp', 'gait_ts_gData', ...
    'gait_ts_eyeposZ_resamp','gait_ts_eyeposY_resamp', ...
    'gait_ts_eyedirZ_resamp', 'gait_ts_eyedirY_resamp', ...
    'doubgait_ts_raw', 'doubgait_ts_resamp','-append');

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Stride (linked gaits?)
 

%     %% sanity check to find long gaits in trials:
%     clf;
%     subplot(211)
%     plot(gait_ts_raw');
%     subplot(212);
%       plot(gait_ts_resamp');
% %%
%     %
%     hold on;
%     for igait = 1:size(gait_ts_raw,1);
%         %plot trial number
%         %plot at end point:
%         tmp = gait_ts_raw(igait,:);
%         lastsamp= max(find(~isnan(tmp)));
%         %if large add text.
%         if lastsamp >75
%             text(lastsamp, tmp(lastsamp), num2str([gait_ts_trialallocation(igait)]))
%         end
%     end % add gait text if long:

end % per ppant

