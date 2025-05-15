%j1A_EyeMovementpertrial.

% Here we will explore eye movements per trial. noting microsaccades and
% any time spent looking off screen (for target/trial rejection).

% plot the trial by trial eye direction data. with interpolated blinks (for
% sanity checks).
%%%%%%

% now updating to use the algorithm implemented in our mobiThings paper.

%% parameters:
visualizeResults=0;     % % useful for debugging

setmydirs_detectv3;

cd(procdatadir)
% show ppant numbers:
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%%
Fs=90;

interpwindow = 6; % iw * 11ms  @ (90 Hz) % this is in smaples.
% interpwindow = .20; % Sec before and after blink to define interpolation.
%%
for ippant = 1:length(pfols)
    cd(procdatadir)


    %     pkdist = participantstepwidths(ippant);
    %%load data from import job.
    load(pfols(ippant).name, 'HeadPos', 'EyePos', 'EyeDir','trialInfo', 'trial_summaryTable', 'subjID');
    savename = pfols(ippant).name;
    disp(['Preparing j1A ' subjID]);

    ftmp = find(savename =='_');
    subjID = savename(1:ftmp(1)-1);

    % try to reorient to fig dir, make subfolder if absent
    try cd([figdir  filesep 'trial_Eye' filesep subjID])
        pfigdir= pwd;
    catch
        mkdir([figdir filesep  'trial_Eye' filesep subjID]);
        cd([figdir filesep  'trial_Eye' filesep subjID])
        pfigdir=pwd;
    end

    %% visualize the peaks and troughs in the head (Y) position.
    %Head position on Y axis, smooth to remove erroneous peaks.
    if visualizeResults
        figure(1); clf;
        set(gcf, 'units', 'normalized', 'position', [0.01,0.01, .9, .9], 'color', 'w', 'visible', 'off');
        %
        pcount=1; % plot indexer
        figcount=1; %figure counter.
    end

    for  itrial=1:size(HeadPos,2)
       
        badtrial=0; % default is that we have a good trial.

        % extract trial data:
        trialtmp = HeadPos(itrial).Y;
        trialD_sm = smooth(trialtmp, 5); % small smoothing factor.
        trial_headData = squeeze(trialtmp);

    [blinksAt, blinksEnd]= deal([]);

        trial_EyeData   = EyeDir(itrial);
        trial_EyeOrigin   = EyePos(itrial);

        Timevec = trialInfo(itrial).times;
        reptrial = find(trial_summaryTable.trial== itrial);
        if isempty(reptrial)
            continue
        end

        isStat = trial_summaryTable.isStationary(reptrial(1));
        isPrac = trial_summaryTable.isPrac(reptrial(1));
        % note some tweaks between participants:
        if iscell(isStat) %convert from cell
            isStat = contains('True', isStat);
        end
        trialTarg = trialInfo(itrial).targstate;
        trialClick = trialInfo(itrial).clickstate;

        %locations of peaks and troughs:
        [locs_ptr, locs_trtr]=deal([]);

        % Note that we can skip when ppant was stationary.

        % detect any blinks.
        % can look for discontinuities as large jumps in the diff.

        %         blinkThreshold varies across participants (default 0.4)

        blinkThresh= adjustBlinkthresh_subjlist_detectver3(subjID);

        Eyetrack = diff(abs([0;trial_EyeData.Y])); % absolute ensures positive diff is blink onset.


        % alternte, find max values.
        EyeX = abs(trial_EyeData.X); % assuming that the blinks will be in both channels.
        EyeY = abs(trial_EyeData.Y);
        EyeZ = abs(trial_EyeData.Z);
%         EyeMax = max([EyeX; EyeY;EyeZ], [],1);

% note that our participants were oriented with Z-Y , walking along the X
% plane. 
        EyeMax = max([EyeZ'; EyeY'], [],1);

        % blinks will be discontinuities (data ==1).
        EyesClosed = find(EyeMax>.80);
       
        if any(EyesClosed)
            % need to convert to samples for the algo below to work:
           
           %automate the interpolation of these clusters.
            %Find clusters:
            A= EyesClosed;
            % Compute the differences between consecutive elements
            diffs = diff(A);

            % Find the indices where the difference is not 1
            breakpoints = find(diffs ~= 1);

            % Compute the start and end indices of each group
            starts = [1, breakpoints+1];
            ends = [breakpoints, numel(A)];
            
            % now repair: if any clusters overlap (within interp window),
            % simply combine.
            blinksAt = EyesClosed(starts);
            blinksEnd= EyesClosed(ends);
            
            % we get problems if the trial started with eyes closed, so
            % remove:
            if blinksAt(1)<interpwindow
                blinksAt=blinksAt(2:end);
                blinksEnd=blinksEnd(2:end);
            end



            
                 % in majority of cases, there is a missing blink end at the
            % trial end. but some subject X trial specific adjustments:
            badtrial=0;
%             adjustBlinkpertrial_mobiThings;
           
            
            
            if badtrial
                subplot(5,3,pcount);
                title('badtrial'); hold on;
                plot(Eyetrack);
                pcount=pcount+1;
                if pcount==16
                    cd(pfigdir)
                    print('-dpng', [subjID ' trialEyeBlinks ' num2str(figcount)  ' ' pwalk]);
                    newFig=1;
                    clf;
                    pcount=1;
                    figcount=figcount+1;
                end
                continue % some trials just not worth the attempt.
            end
          

            blinksDur = blinksEnd- blinksAt;
            % %define our interpretation window (with a buffer)
%             bufferFrom= dsearchn(Timevec, [Timevec(blinksAt)- 2*interpwindow]');
% 
%             bufferUntil= dsearchn(Timevec, [Timevec(blinksEnd)+ 2*interpwindow]');

            bufferFrom= blinksAt- 2*interpwindow;

            bufferUntil= blinksEnd+ 2*interpwindow;
    % avoid under/overshoot
            bufferFrom(bufferFrom<=0)=1;
            bufferUntil(bufferUntil>length(trialTarg))=length(trialTarg);
            
            %% interpolate
            for iblink = 1:length(bufferFrom)
                datachunk_Y = trial_EyeData.Y(bufferFrom(iblink):bufferUntil(iblink));
                datachunk_Z = trial_EyeData.Z(bufferFrom(iblink):bufferUntil(iblink));
                datachunk_X = trial_EyeData.X(bufferFrom(iblink):bufferUntil(iblink));
                % also interp eye origin info:
                datachunk_oX= trial_EyeOrigin.X(bufferFrom(iblink):bufferUntil(iblink));
                datachunk_oY= trial_EyeOrigin.Y(bufferFrom(iblink):bufferUntil(iblink));
                datachunk_oZ= trial_EyeOrigin.Z(bufferFrom(iblink):bufferUntil(iblink));
    
                %now determine the window to interp, within this buffered
                %epoch:
                % start  interp window before the blink in this
                % epoch:
                adjustedOnset = interpwindow; % WE HAD A 2*interp window buffer before blink start.
                if adjustedOnset==0
                    adjustedOnset=1;
                end
                % end interp window after interpwindow at the blink end in
                % this epoch.

                adjustedOffset = length(datachunk_Y) - interpwindow; 
                blinkwindow = adjustedOnset:adjustedOffset;

                %%
                % interpolate all data
                %  YQ = spline(X,Y,XQ) performs cubic spline interpolation using the
                %     values Y at sample points X to find interpolated values YQ at the query
                %     points XQ.
                %         - X must be a vector.
                %         - If Y is a vector, Y(j) is the value at X(j).
                %         - If Y is a matrix or n-D array, Y(:,...,:,j) is the value at X(j).
                %
                x= 1:length(datachunk_X);
                xq = blinkwindow;
                repairData= {datachunk_X, datachunk_Y, datachunk_Z, datachunk_oX, datachunk_oY,datachunk_oZ};
                cleanData=[];
                for id= 1:length(repairData)
                    reptmp = repairData{id};
                    reptmp(blinkwindow) =nan;

                    % these methods may introduce overshoot:
%                     YQ = spline(x,reptmp,xq);
%                     YQ = makima(x, reptmp, xq);

%safest probably just to linearly interpolate:

                    % Find indices of missing data
                    nan_indx = isnan(reptmp);
%                     %replace nan with 0 for interp1 below.
                    v_interp= reptmp;
                    v_interp(nan_indx)=0;
% 
%                     % alternate Interpolate missing data
                    YQ = interp1(find(~nan_indx), v_interp(~nan_indx), 1:numel(v_interp), 'linear'); % use linkear polynomial.

                    % place new data in original
%                     datachunk_clean = repairData{id};
%                     datachunk_clean(blinkwindow)= YQ;

                    datachunk_clean=YQ;
                    cleanData{id}= datachunk_clean; % for plots

                    %now we have the repaired data, replace in EyeDir:
                    switch id
                        case 1
                            %now we have the repaired data, replace in EyeDir:
                            trial_EyeData.X(bufferFrom(iblink):bufferUntil(iblink))= datachunk_clean;
                        case 2
                            trial_EyeData.Y(bufferFrom(iblink):bufferUntil(iblink))= datachunk_clean;
                        case 3
                            trial_EyeData.Z(bufferFrom(iblink):bufferUntil(iblink))= datachunk_clean;

                        case 4
                            trial_EyeOrigin.X(bufferFrom(iblink):bufferUntil(iblink))= datachunk_clean;
                        case 5
                            trial_EyeOrigin.Y(bufferFrom(iblink):bufferUntil(iblink))= datachunk_clean;
                        case 6
                              trial_EyeOrigin.Z(bufferFrom(iblink):bufferUntil(iblink))= datachunk_clean;

                    end
                end

                %% sanity check the spline worked:
% %                                 figure(10); clf
%     titlesare={'x', 'y', 'z', 'Xorigin'};
%                     for id=1:4
%     
%                     subplot(2,2,id)
%                     plot(repairData{id}, 'ko-', 'linew', 1); hold on;  shg
%                     for ip=1:length(blinkwindow)                        
%                     plot(blinkwindow(ip), repairData{id}(blinkwindow(ip)), 'ro')
%                     end
%                     title(titlesare{id});
%                 end
            end % per blink in trial:


        end % any eye closed.


       
        %we may want to vis the trial onsets and response:
        if visualizeResults

            % small function to plot per panel:
            plotD=[];
            plotD.itrial=itrial;
            plotD.figcount=figcount;
            plotD.Eye_O_clean= trial_EyeOrigin;
            plotD.Eye_D_clean= trial_EyeData;
            plotD.Eye_O_raw= HeadPos(itrial);
            plotD.Eye_D_raw= EyeDir(itrial);
            %show where interp happened:
            plotD.blinksAt=blinksAt;
            plotD.blinksEnd=blinksEnd;
             plotD.interpWindow= interpwindow; % show interp boundaries.

            plotD.Timevec= Timevec;
            
            plotD.pfigdir= pfigdir;
            plotD.subjID= subjID;

%             plotD.walkSpeed= pwalk;
%              plotD.walkSpeedInt= wlkSpeed;

            % see funtion defn below.
            newFig=quickplot(pcount,plotD);
           
            
            if newFig==1
                pcount=2;
                figcount=figcount+1;
            else
                pcount=pcount+1;
                if pcount==16

                    cd(pfigdir)
                    print('-dpng', [subjID ' trialEyeBlinks ' num2str(figcount)  ]);
                    newFig=1;
                    clf;
                    figcount=figcount+1;
                    pcount=1;
                end

            end
        end % if visualizing...

        %% store results.
        EyeDir(itrial).Y_clean = trial_EyeData.Y;
        EyeDir(itrial).Z_clean = trial_EyeData.Z;
        EyeDir(itrial).X_clean = trial_EyeData.X;
         EyeDir(itrial).blinksAt = blinksAt;
        EyeDir(itrial).blinksEnd= blinksEnd;
        

        EyePos(itrial).Y_clean = trial_EyeOrigin.Y;
        EyePos(itrial).Z_clean = trial_EyeOrigin.Z; % we will call later to plot gait_cycle effects.
        EyePos(itrial).X_clean = trial_EyeOrigin.X; 
        EyePos(itrial).blinksAt = blinksAt;
        EyePos(itrial).blinksEnd= blinksEnd;

        HeadPos(itrial).badtrial= badtrial; % hard skip some cases.


    % note we can also calculate angular velocity for saccade detection:
    % based on algo of Nystrom and holmqvist:

%% plots the euclidean velocity:
% clf
% tmp=diff(trial_EyeData.Y);
% plot(tmp, 'k');
% fltord = 2;
% fltlen= 5; % 2 x min saccade duration (min = 200 ms 2frames)
% B= sgolayfilt(tmp,fltord, fltlen);
% hold on;
% plot(B, 'r', 'LineWidth',2)
% C= sgolayfilt(diff(trial_EyeData.Z),fltord, fltlen);
% plot(C,'b');
% % euclidean distance between points:
% eD=[];
% for isamp= 1:length(B);
% eD(isamp) = norm(B(isamp) - C(isamp));
% end
% clf
% plot(eD, 'm');
% hold on;
% overlay targets to see if there is a saccade?


%%
    end %itrial.
    %% after all trials, print the final figure (in case uneven subplots/trial counts).

   if visualizeResults
    cd(pfigdir)
    print('-dpng', [subjID ' trialEyeBlinks ' num2str(figcount) ]);
end
    %resave with new data in structure.
    cd(procdatadir)
    save(savename, 'EyeDir', 'EyePos', '-append');

end % isub
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calling this function above.
function newfig= quickplot(pcount, plotD)
% function to handle subplot placement, and printing/counting new figures.
%newfig is a flag to start the next figure, toggles subplot count.
figure(1);
newfig=0;
if pcount ==16
    %print that figure, and reset trial count.

    cd(plotD.pfigdir)
    print('-dpng', [plotD.subjID ' trialEyeBlinks ' num2str(plotD.figcount)  ' ' plotD.walkSpeed]);
    newfig=1;
    clf;
    pcount=1;
end
subplot(5,3,pcount);

plot(plotD.Timevec, plotD.Eye_D_raw.Z, 'k'); hold on;
plot(plotD.Timevec, plotD.Eye_D_raw.Y, 'k');
% plot(plotD.Timevec, plotD.Eye_D_raw.X, 'r');
xh=plot(plotD.Timevec, plotD.Eye_D_clean.Z, 'b-'); hold on;
yh=plot(plotD.Timevec, plotD.Eye_D_clean.Y, 'b-');
axis tight;


% plot(plotD.Timevec, plotD.Eye_D.X, 'r');

ylabel('Eye direction (position)');
% ylim([-.5 .5]);
iw= plotD.interpWindow;
for iblink = 1:length(plotD.blinksAt);
%     plot([plotD.blinksAt(iblink), plotD.blinksAt(iblink)], ylim, 'k', 'LineWidth',2)
%         plot([plotD.blinksEnd(iblink), plotD.blinksEnd(iblink)], ylim, 'k')
xvec=[plotD.blinksAt(iblink), plotD.blinksAt(iblink), plotD.blinksEnd(iblink), plotD.blinksEnd(iblink)];
xvec = plotD.Timevec(xvec);
yvec=[-.5 .5 .5 -.5];

ph=patch(xvec,yvec, [.2 .2 .2]);
ph.FaceAlpha= .2;
ph.EdgeColor= 'w';
% pg.
end




yyaxis right;
plot(plotD.Timevec, plotD.Eye_O_clean.Z - mean(plotD.Eye_O_clean.Z), 'r:'  ); hold on;
plot(plotD.Timevec, plotD.Eye_O_clean.Y- mean(plotD.Eye_O_clean.Y) , 'b:'  ); % normalized to account for different heights.
hold on;
ylabel('Eye origin (m)');
xlabel('Time (s)');
title(['Trial ' num2str(plotD.itrial) ]);
set(gca,'YColor', 'k')
legend([xh,yh], {'X(Z)','Y'})
% disp(['STD of trial ' num2str(plotD.itrial) ' ' sprintf('%.3f',std(plotD.Eye_D.Z))])
% axis tight
end

