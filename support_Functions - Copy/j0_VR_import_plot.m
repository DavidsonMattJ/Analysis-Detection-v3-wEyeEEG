% Walking + VR contrast detection experiment v3.
% This script imports the raw csv files(*FramebyFrame, *trial_summary) and 
% resaves after some basic preprocessing and wrangling for further Matlab 
% analysis.
%% set paths and directories.

setmydirs_detectv3;
cd(datadir)
sancheckplots=0; % toggle to show sanity check plots to check input.
pfols = dir([pwd filesep '*framebyframe.csv']);
nsubs= length(pfols);

% show ppant numbers: in command window
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%% Per csv file, import and wrangle into Matlab Structures, and data matrices:
for ippant =1:length(pfols)
    cd(datadir)

    pfols = dir([pwd filesep '*framebyframe.csv']);

    %% load subject data as table.
    filename = pfols(ippant).name;
    %extract name&date from filename:
    ftmp = find(filename =='_');
    subjID = filename(1:ftmp(end)-1);


    savename = [subjID '_summary_data'];

    %query whether we want to recompute frame x frame (unlikely).
    cd(procdatadir);
    %     if exist([savename '.mat'], 'file')==2
    %         disp(['frame x frame alread saved for ' subjID]);
    %     else
    %read table
    cd(datadir);
    opts = detectImportOptions(filename,'NumHeaderLines',0);
    disp(['reading large frame x frame file now...']);
    T = readtable(filename,opts);
    ppant = T.participant{1};
    disp(['Preparing participant ' ppant]);

    [TargPos, HeadPos, TargState, ClickState, EyePos, EyeDir] = deal([]);

    %% use logical indexing to find all relevant info (in cells)
    posData = T.position;
    clickData = T.clickstate;
    targStateData= T.targState;

    objs = T.trackedObject;
    axes= T.axis;
    Trials =T.trial;
    Times = T.t;

    targ_rows = find(contains(objs, 'target'));
    head_rows = find(contains(objs, 'head'));
    eyePos_rows = find(contains(objs, 'gazeOrigin'));
    eyeDir_rows = find(contains(objs, 'gazeDirection'));

    Xpos = find(contains(axes, 'x'));
    Ypos  = find(contains(axes, 'y'));
    Zpos = find(contains(axes, 'z'));
    %%
    userows = {head_rows, targ_rows, eyePos_rows, eyeDir_rows};
    for idatatype = 1:length(userows)

        %% per type, find the intersect of thse indices, to fill the data.
        datarows = userows{idatatype};
        Dx = intersect(datarows, Xpos);
        Dy = intersect(datarows, Ypos);
        Dz = intersect(datarows, Zpos);

        %% further store by trials (walking laps).
        vec_lengths=[];

        DataPos=[]; % will be renamed below.

        for itrial = 1:length(unique(Trials))

            trial_rows = find(Trials==itrial-1); % indexing from 0 in Unity

            DataPos(itrial).X = posData(intersect(Dx, trial_rows));
            DataPos(itrial).Y = posData(intersect(Dy, trial_rows));
            DataPos(itrial).Z = posData(intersect(Dz, trial_rows));

            % only need to perform once, but also capture the targstate and
            % click state, and times on each trial

            trial_times = Times(intersect(Dx, trial_rows));

            if idatatype==1
                trialInfo(itrial).targstate= targStateData(intersect(Dx, trial_rows));
                trialInfo(itrial).clickstate= clickData(intersect(Dx, trial_rows));
                trialInfo(itrial).times = trial_times;

                %                 plot(trialInfo(itrial).clickstate); hold on;

            end
        end

        % save per datatype:
        switch idatatype
            case 1
                HeadPos = DataPos;
            case 2
                TargPos = DataPos;
            case 3
                EyePos = DataPos;
            case 4
                EyeDir= DataPos;

        end
    end % idatatype
    %%
    disp(['Saving position data split by trials... ' subjID]);
    cd(procdatadir)
    try save(savename, 'TargPos', 'HeadPos', 'EyePos', 'EyeDir', 'trialInfo', 'subjID', 'ppant', '-append');
    catch
        save(savename, 'TargPos', 'HeadPos', 'EyePos', 'EyeDir','trialInfo', 'subjID', 'ppant');
    end

    %     end
    %
    %% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ now summary data
    cd(datadir)
    pfols = dir([pwd filesep '*trialsummary.csv']);
    nsubs= length(pfols);
    %
    filename = pfols(ippant).name;

    %extract name&date from filename:
    ftmp = find(filename =='_');
    subjID = filename(1:ftmp(end)-1);
    %read table
    opts = detectImportOptions(filename,'NumHeaderLines',0);
    T = readtable(filename,opts);
    rawSummary_table = T;
    disp(['Preparing participant ' T.participant{1} ]);


    savename = [subjID '_summary_data'];
    % summarise relevant data:
    targPrestrials =find(T.nTarg>0);
    practIndex = find(T.isPrac ==1);
    npracTrials = (T.trial(practIndex(end)) +1);
    disp([subjID ' has ' num2str(npracTrials) ' practice trials']);

    %extract the rows in our table, with relevant data for assessing
    %calibration
    nstairs = unique(T.qStaircase(T.qStaircase >0));
    calibAcc=[];
    for iqstair=1:length(nstairs) % check all staircases.
        tmpstair = nstairs(iqstair);
        qstairtrials= find(T.qStaircase==tmpstair);
        calibAxis = intersect(qstairtrials, practIndex);

        %calculate accuracy (running accuracy).
        calibData = T.targCor(calibAxis);
        for itarg=1:length(calibData)
            tmpD = calibData(1:itarg);
            calibAcc(iqstair).d(itarg) = sum(tmpD)/length(tmpD);
        end

    end

    %% FAlse alarms: repair FA data in table (happens if we have multiple FAs, extra columns are added, which need to be collapsed.
    ColIndex = find(strcmp(T.Properties.VariableNames, 'FA_rt'), 1);
    %repair string columns:
    for ix=ColIndex:size(T,2)
        tmp = table2array(T(:, ix));
        if iscell(tmp)
            % where is there non-zero data, to convert?
            replaceRows= find(~cellfun(@isempty,tmp));
            %convert each row
            newD=zeros(size(tmp,1),1);
            for ir=1:length(replaceRows)
                user = replaceRows(ir);
                newD(user) = str2num(cell2mat(tmp(user )));
            end

            %% change table column type:
            colnam = T.Properties.VariableNames(ix);
            T.(colnam{:}) = newD;
            %             T(:,ix) =table(newD');

        end


    end
    %% extract Target onsets per trial (as table).
    %% and Targ RTs, contrast values if they exist.
    alltrials = unique(T.trial);

    %we want to repair the table as we go.
    trial_summaryTable= T;
    %remove some cols we don't need:
    trial_summaryTable.date=[];
    trial_summaryTable.nTarg=[];
    trial_summaryTable.qStaircase=[];
    trial_summaryTable.intActualE=[];

    %%
    for itrial= 1:length(alltrials)
        thistrial= alltrials(itrial);
        relvrows = find(T.trial ==thistrial); % unity idx at zero.

        %create vectors for storage:
        tOnsets = T.targOnset(relvrows);
        tRTs = T.targRT(relvrows);
        tCor = T.targCor(relvrows);
        tFAs= table2array(T(relvrows(end), ColIndex:end));
        tFAs= tFAs(~isnan(tFAs));

        %% Reaction time region (tidy / reclassify some values).

        RTs = (tRTs - tOnsets);
        %         %note that negative RTs, indicate that no response was recorded:
        tOmit = find(RTs<=0);
        if ~isempty(tOmit)
            %             tCor(tOmit) = NaN; % don't count those incorrects, as a mis identification.
            RTs(tOmit)=NaN; % remove no respnse
        end

        % now we can add this RT information to our revised table
        trial_summaryTable.clickRT(relvrows) = RTs;


        %% we also want to reject very short RTs (reclassify as a FA).

        if any(find(RTs<0.15))
            disp(['Trial: ' num2str(thistrial) ' Suspicious RT']);
            checkRT= find(RTs<.15);
            %debug to check:

            % reclassify data (mark as incorrect, and record as FA).
            for ifalseRT = 1:length(checkRT)
                indx = relvrows(checkRT(ifalseRT));
                trial_summaryTable.targCor(indx)= 0;
                trial_summaryTable.clickRT(indx)= NaN;
                trial_summaryTable.FA_rt(indx) =  trial_summaryTable.targRT(indx);

                %debug add to fig:
                if sancheckplots
                    clf; title(['Trial: ' num2str(thistrial) 'early RT?']); hold on;
                    plot(trialInfo(itrial).times, trialInfo(itrial).targstate); hold on;
                    plot(trialInfo(itrial).times, trialInfo(itrial).clickstate);
                    %
                    plot([tRTs(checkRT(ifalseRT)), tRTs(checkRT(ifalseRT))], [0 1], 'om')
                end
            end
            %             pause; % wait for click

        end
        %%
        % seems some FA are missing, do a quick check to see if
        % any extra in the pure clickdata.
        tRTs= tRTs(tRTs>0);
        clks= find(trialInfo(itrial).clickstate);
        clksTs= trialInfo(itrial).times(clks);
        if length(clksTs) ~= length(tRTs)
            %FA present
            %find the outlier
            %%
            %sanity debug:
            if sancheckplots
                clf;
                title(['Trial: ' num2str(thistrial) ' False Alarms? (frXfr data)']); hold on;
                tg=plot(trialInfo(itrial).times, trialInfo(itrial).targstate, 'b'); hold on;
                cl=plot(trialInfo(itrial).times, trialInfo(itrial).clickstate, 'r');


                hold on
                %also plot the target onsets and recorded RTs, to see which are missing:
                for irt = 1:length(tOnsets)
                    tgl = plot([tOnsets(irt) tOnsets(irt)], ylim, 'o-b');
                end
                for irt = 1:length(tRTs)
                    rtl=plot([tRTs(irt) tRTs(irt)], ylim, 'o-r');
                end
                try
                    legend([tg, cl, tgl, rtl], {'targState', 'clickState', 'tsumry', 'cksmry'})
                catch
                end
                shg
            end
            %%
            %round times to 3dp:
            clksTsR= round(clksTs,3);
            % find out which is furthest from a recorded click
            % in the sumry data. .:. a FA:
            if ~isempty(tRTs)
                [loc, dist] =dsearchn(tRTs, clksTsR);
                FAidx= find(dist > .1);
            else
                FAidx = 1:length(clksTsR);
            end
            for iextraclick = 1:length(FAidx)
                thisFA_isat = clksTsR(FAidx(iextraclick));

                % note, if this FA isn't already recorded,
                % store
                if ~isempty(tFAs)
                    [~, distfromrecd] = dsearchn(tFAs', thisFA_isat);
                    if distfromrecd >.05 % another click
                        tFAs= [tFAs, thisFA_isat];
                    end

                else % first FA
                    tFAs =thisFA_isat;
                end
            end

            %avoid junk FA and first frame.
            tFAs= tFAs(tFAs>0);

            % add our new FAs to the revised table:
            nantmp = nan(1, length(relvrows));
            nantmp(1:length(tFAs)) =tFAs;
            if length(relvrows)~= length(nantmp)
                %. % omits last FAs, but shouldn't be a problem, since we
                %don't use last gaits in each trial.
                trial_summaryTable.FA_rt(relvrows) = nantmp(1:length(relvrows));
            else
                trial_summaryTable.FA_rt(relvrows) = nantmp;
            end
            %add FA to fig:
            if sancheckplots
                for irt = 1:length(tFAs)
                    rtl=plot([tFAs(irt) tFAs(irt)], ylim, 'o-m');
                end

                ylim([0 1.5])
            end
            %              pause; % await click
        end


        % end RT region


    end
    %repair trial indexing:
    trial_summaryTable.trial =  trial_summaryTable.trial +1;
    trial_summaryTable.block =  trial_summaryTable.block +1;
    trial_summaryTable.trialID =  trial_summaryTable.trialID +1;

    %also convert cells to easier format (only some ppants).
    isstatall= trial_summaryTable.isStationary;

    if iscell(isstatall)

        isstatint = contains(isstatall, 'TRUE', 'IgnoreCase',true);
        trial_summaryTable.isStationary = isstatint;
    end

    %save this structure for later analysis per gait-cycle:
    disp(['Saving trial summary data ... ' subjID]);
    cd(procdatadir)
    save(savename, 'trial_summaryTable','-append');


end % participant
