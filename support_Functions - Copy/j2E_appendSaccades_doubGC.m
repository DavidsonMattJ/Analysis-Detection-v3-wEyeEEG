% j2E_appendSaccades_doubGC

% building off the back of j2D- > adding an extra column for the
% classification of an event in either a LRL or RLR step sequence.
% necessary for showing double step proportions of DVs etc.

% this reclassifies saccade onsets as a proportion of the stride-cycle.

% quest detect v3 ! 

%laptop:
setmydirs_detectv3;

cd(procdatadir)
%% show ppant numbers:
pfols = dir([pwd filesep '*_summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%
%%


for ippant= 1:nsubs
    cd(procdatadir)    %%load data from import job.   
    load([pfols(ippant).name], 'saccadeSummaryTable', 'HeadPos');
    savename = pfols(ippant).name;
    disp(['Preparing j2 cycle data... ' savename]);

    
    % we will need to call the above head pos based on event info:

    %% Gait extraction.
    % per event, find the relevant gait. then add the gPcnt relative
    % to LRL anr RLR sequence (if possible).

    allevents = size(saccadeSummaryTable,1);

    for ievent= 1:allevents

        %for each event (row), convert to pcnt of linked steps (stride
        %cycle).

        tTrial = saccadeSummaryTable.trialID(ievent);

        
        trs = HeadPos(tTrial).Y_gait_troughs;

        thisOnset = saccadeSummaryTable.saccadeOnsetFrame(ievent);
        thisOffset = thisOnset + saccadeSummaryTable.saccDur_samps(ievent);

        %for onset and offset, find percent relative to double gc.
        % find nearest gait
        gindx_Onset = saccadeSummaryTable.gindx_Onset(ievent);
        gindx_Offset = saccadeSummaryTable.gindx_Offset(ievent);


        % extract the event as % gait.[0 100]
        useINDX = [gindx_Onset, gindx_Offset];
        useFT =  {saccadeSummaryTable.saccOns_gFoot(ievent),saccadeSummaryTable.saccOffs_gFoot(ievent)};
        useSamp= [thisOnset,thisOffset];
        saveCols= {'saccOnset', 'sacOffset'};
        for igOnOff=1:2
            % sacconset:
            gindx = useINDX(igOnOff);
          
            % we have the gIndx, so look at prev and current step.
            %
            currentFt = useFT{igOnOff}; % note this was provided in earlier job 2D
          
            thisOnset= useSamp(igOnOff);
            if strcmp(currentFt, 'LR')
                currDft = 'RLR';
                prevDft = 'LRL';
            else
                currDft = 'LRL';
                prevDft= 'RLR';
            end
            
            % 
            gPcnt_step1=nan;
            gPcnt_step2=nan;
            
             % % extract the event as % double gait.[0 100]
            % might have been over complicating things.
            % we have the current gait index. there are 2 strides to account for. 
            % this onset within the first step of a stride, and same onset
            % within the second step of a stride (i.e. 2 percentages to
            % calculate).
            
            % guard clause for the end cases (e.g. before stride has begun).
            

            %first prcntage is easy, simply stretch to a proportion of the
            %time for current trough to trough+2
            if gindx<=(length(trs)-2)
            firstSamps= trs(gindx):trs(gindx+2);
%             % percentage of this onset as a proportion of this interval:
             tDur = firstSamps(end)- firstSamps(1);
             tE= thisOnset-firstSamps(1);
             gPcnt= ceil((tE/tDur)*100); % as a prcntage of the whole stride:
             gPcnt_step1 = gPcnt;
%                 
%              %else original way (single step
%                firstSamps= trs(gindx):trs(gindx+1);
%                 tDur = firstSamps(end)- firstSamps(1);
%              tE= thisOnset-firstSamps(1);
%              gPcnt_step1= round((tE/tDur)*50);
%              

            elseif gindx==(length(trs)-1) % perform estimation:
                % we can't go forward 2 steps, so we can estimate:

                firstSamps= trs(gindx):trs(gindx+1);
            % percentage of this onset as a proportion of this interval:
             tDur = firstSamps(end)- firstSamps(1);
             tE= thisOnset-firstSamps(1);
             gPcnt_step1= round((tE/tDur)*50); % note the denominator.


            end

            if gPcnt_step1>50 
%                 error('check step1')
            end
%              now calculate as a proportion of the second step in a
%              stride:
            if gindx>1  &&  gindx<=length(trs)-1% otherwise we can skip/do other.
             secondSamps= trs(gindx-1):trs(gindx+1);
              tDur = secondSamps(end)- secondSamps(1);
             tE= thisOnset-secondSamps(1);
             gPcnt= ceil((tE/tDur)*100);
             gPcnt_step2 = gPcnt;


% %og way.
%                 secondSamps= trs(gindx):trs(gindx+1);
%               tDur = secondSamps(end)- secondSamps(1);
%              tE= thisOnset-secondSamps(1);
%              gPcnt= round((tE/tDur)*50);
%              gPcnt_step2 = gPcnt;
%              gPcnt_step2 = gPcnt+50;

            elseif gindx==1 % need to estimate the position as a second step location.
                % use the duration of the current step.
                secondSamps = trs(gindx):trs(gindx+1);
                tDur = secondSamps(end)- secondSamps(1);
                tE= thisOnset-secondSamps(1);
                gPcnt= round((tE/tDur)*50);
                %
                % rather than using 50 outright, 
                % use the mean step length on this trial:

%                 gPcnt_step2 = gPcnt+50;

                 gPcnt_step2 = gPcnt+50;%round(median(diff(trs)));

            end

            if gPcnt_step2>100
                gPcnt_step2=nan;
                disp('check step 2')
            end

             %store .
            
%            
            %>>>>> add this new info to our table,
            saccadeSummaryTable.([saveCols{igOnOff}  '_gPcnt_' currDft])(ievent)= gPcnt_step1;
            saccadeSummaryTable.([saveCols{igOnOff}  '_gPcnt_' prevDft])(ievent)= gPcnt_step2;
            saccadeSummaryTable.([saveCols{igOnOff}  '_gPcnt_step1inStride'])(ievent)= gPcnt_step1;
            saccadeSummaryTable.([saveCols{igOnOff}  '_gPcnt_step2inStride'])(ievent)= gPcnt_step2;


        end % onset offset




    end % each row in table (event)

    % debug. Show trgOs
    % ptrgOs = trial_summaryTable.trgO_gPcnt;
    % ptrgOs = ptrgOs(~isnan(ptrgOs));
    % hist(ptrgOs, 100);

  
    %%
    %%
    disp(['Finished J2E appending DOUBLE gait percentage data for ... ' subjID]);
    save( [pfols(ippant).name], 'saccadeSummaryTable', '-append'); % overwrite.
end % participant

