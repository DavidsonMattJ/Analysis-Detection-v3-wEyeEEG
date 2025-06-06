% AA_checkEyeData\

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
% wrangle frame by frame.
for ippant=1
% cd(datadir)
    
    pfols = dir([pwd filesep '*framebyframe.csv']);

    %% load subject data as table.
    filename = pfols(ippant).name;
    %extract name&date from filename:
    ftmp = find(filename =='_');
    subjID = filename(1:ftmp(end)-1);
    
   
    savename = [subjID '_summary_data'];
    
    %query whether pos data job has been done (is in list of variables
    %saved)
%     cd('ProcessedData')
%     listOfVariables = who('-file', [savename '.mat']);
%     if ~ismember('HeadPos', listOfVariables)  
        % if not done, load and save frame x frame data.
    % simple extract of positions over time.
    
    %read table
    opts = detectImportOptions(filename,'NumHeaderLines',0);
    T = readtable(filename,opts);
    ppant = T.participant{1};
    disp(['Preparing participant ' ppant]);
    
    [TargPos, HeadPos,  EyePos, EyeDir,TargState, clickState] = deal([]);
    
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
   
    eye_rows = find(contains(objs, 'gazeOrigin'));
    eyedir_rows = find(contains(objs, 'gazeDirection'));
    
    Xpos = find(contains(axes, 'x'));
    Ypos  = find(contains(axes, 'y'));
    Zpos = find(contains(axes, 'z'));
    
    objs = {targ_rows, head_rows, eye_rows, eyedir_rows};
    
    for iobj=1:4
        useobjrows= objs{iobj};
    %% now find the intersect of thse indices, to fill the data.
    hx = intersect(useobjrows, Xpos);
    hy = intersect(useobjrows, Ypos);
    hz = intersect(useobjrows, Zpos);
    
    
    %% further store by trials (walking laps).
    vec_lengths=[];
    objPos=[];
    for itrial = 1:length(unique(Trials))
        
        trial_rows = find(Trials==itrial-1); % indexing from 0 in Unity
        
        trial_times = Times(intersect(hx, trial_rows));
        %Head first (X Y Z)
        objPos(itrial).X = posData(intersect(hx, trial_rows));
        objPos(itrial).Y = posData(intersect(hy, trial_rows));
        objPos(itrial).Z = posData(intersect(hz, trial_rows));
        %store time (sec) for matching with summary data:
        objPos(itrial).times = trial_times;
        objPos(itrial).isStationary = unique(T.isStationary(trial_rows));
        
        
        % because the XYZ have the same time stamp, collect click and targ
        % state as well.
        % note varying lengths some trials, so store in structure:
        if iobj==1
        TargState(itrial).state = targStateData(intersect(hx, trial_rows));
        TargState(itrial).times = trial_times;
        clickState(itrial).state = clickData(intersect(hx, trial_rows));
        clickState(itrial).times = trial_times;
        end
        
        
        %switch and rename

        
    end
    switch iobj %TargPos, HeadPos,  EyePos, EyeDir,
        case 1
            TargPos= objPos;
        case 2
            HeadPos= objPos;
        case 3
            EyePos = objPos;
        case 4
            EyeDir = objPos;
    end
    end
%% plot walk:





fontSize=10;
for itrial = 18
    
    trialTargPos = TargPos(itrial);
    trialEyePos = EyePos(itrial);    
    trialEyeDir = EyeDir(itrial);
    trialHeadPos = HeadPos(itrial);
   plotD={trialTargPos, trialEyePos, trialEyeDir, trialHeadPos};
   titles={'TargPos', 'EyePos', 'EyeDir', 'HeadPos'};
   
    figure(1); clf
    ic=1;
    for id=1:4
    useD= plotD{id};
    subplot(4,4,ic);ic=ic+1; hold on;
    plot(useD.X, 'k');    
    plot(useD.Y, 'b');    
    plot(useD.Z, 'r'); title([ titles{id} ' XYZ overlayed']);
    
    subplot(4,4,ic); ic=ic+1;
    plot(useD.X, 'k');    title('X'); 

    subplot(4,4,ic); ic=ic+1;
    plot(useD.Y, 'b');    title('Y'); 
    subplot(4,4,ic);ic=ic+1;
    plot(useD.Z, 'r');     title('Z');
    end
end
%%     %%
numberOfSteps = length(trialTargPos.X);

xy = zeros(numberOfSteps,2);
%
clf
plottimes =  [0,1,2,3,4,5,6,7,8];
texttimes = [2; dsearchn([TargState(itrial).times],  [1,2,3,4,5,6,7,8]')];

posData = {trialTargPos, trialEyePos, trialEyeDir};
useCols = {'b', 'k', 'm'};
for postoPlot=[1,3]
    subplot(3,1,postoPlot);
    pdata = posData{postoPlot};
ic=1;
for iframe = 2 : numberOfSteps
	% Walk in the x direction.
	
	% Now plot the walk so far.
	xCoords = pdata.Z(iframe-1:iframe);
	yCoords = pdata.Y(iframe-1:iframe);
    if TargState(itrial).state(iframe)==1
        
        linespecs = [useCols{postoPlot} 'o-'];
        LineWidth = 3;
        
    else
        linespecs = [useCols{postoPlot} '.'];        
        LineWidth = 1;
    end
        
	plot(xCoords, yCoords,  linespecs, 'LineWidth', LineWidth);
	hold on;
    
    if ismember(iframe, texttimes)
        text(xCoords,yCoords, [num2str(plottimes(ic))]);
    ic=ic+1;    
    end
end %each frame
axis tight
end %datatype
%%

end %itrial


% end % ippant


