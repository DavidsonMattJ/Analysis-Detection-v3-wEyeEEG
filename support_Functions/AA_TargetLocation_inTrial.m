% AA_TargetLocation_inTrial

%%%%%% QUEST & EEG version

% Detection experiment (contrast)
%%  Import from csv. FramebyFrame, then summary data.

%laptop:
setmydirs_detectv3
cd(procdatadir);
pfols= dir([pwd filesep '*_summary_data.mat']);
%%
%various jobs as sanity checks


jobs=[]; % 
jobs.concat_plotTargetLocation =1; % all target locations, at onset.

GFX_allTargLocs=[];
for ippant=1:length(pfols)
    % cd(datadir)

    pfols = dir([pwd filesep '*summary_data.mat']);

    %% load subject data as table.
    filename = pfols(ippant).name;
    %extract name&date from filename:
    ftmp = find(filename =='_');
    subjID = filename(1:ftmp(end)-1);


    savename = [subjID '_summary_data'];

    load(filename, 'trial_summaryTable', 'TargPos', 'trialInfo', 'EyeDir', 'EyePos');

    %% plot walk:

    % for each walk type (7 or 8). plot distribution of trial target onsets.
    alltrials = unique(trial_summaryTable.trial);
    allX=[];
    allY=[];
    eyeX=[];
    eyeY=[];
    allOffset_theta=[];
    allOffset_DVA=[];
    allOffset_DVA2=[];

    allFlag_rej=[];
    
    for itrial = 1:length(alltrials)

        usetrial= alltrials(itrial);
        user= find(trial_summaryTable.trial == usetrial);
        tOns = trial_summaryTable.targOnset(user);
        % adapt to check DVA at targ offset (so after movement to target).
        tOns= tOns+.02; % 20 ms duration

        ttimes = dsearchn(trialInfo(usetrial).times, tOns);
        tPosX = round(TargPos(usetrial).Z(ttimes),3);
        tPosY = round(TargPos(usetrial).Y(ttimes),3);
        % and eye:
        eyePosX = round(EyeDir(usetrial).Z_clean(ttimes),3);
        eyePosY = round(EyeDir(usetrial).Y_clean(ttimes),3);

%subtract mean to centre at 0.
tPosX= tPosX- mean(TargPos(usetrial).Z); % of all points.
tPosY= tPosY- mean(TargPos(usetrial).Y); % of all points.



        allX= [allX; tPosX];        
        allY= [allY; tPosY];


        eyeX = [eyeX; eyePosX];
        eyeY = [eyeY; eyePosY];



    % using the HMD distance from target, we can calculate the offset as
    % dva.
    hO = EyePos(usetrial).X(ttimes);
    hDist = abs(EyePos(usetrial).X(ttimes) - TargPos(usetrial).X(ttimes));

%hDist is in metres.
    for itarg= 1:length(tPosX)

%may also be interested in the angle. 

% convert target to origin:
X= [eyePosX(itarg)-tPosX(itarg)];
Y= [eyePosY(itarg)-tPosY(itarg)];
theta = atan2(Y,X);
rho = sqrt(X^2 + Y^2);

% using the distance, convert to DVA.
% pixel relative to distance and width of screen.
% pixelSize = hDist(itarg)*tan(pi/180) / W;
% gaze offset from target in dva:
% offset_degrees = atan(rho * pixelSize / hDist(itarg)) * 180 / pi;
 %% alternatively: 
 % calculate the vector as a dva:

sz = rho; % size of the object (here the offset, in m).

ang = 2*180*atan(sz/(2*hDist(itarg)))/pi;

sz_Ann = .11; % fixed at 11 cm diameter
angAnn =  2*180*atan(2*sz_Ann/(2*hDist(itarg)))/pi;

% flag for rejection if offset exceeds diameter of anulus.
tmpRej = ang>angAnn;
 %%

allOffset_theta=[allOffset_theta, theta];
% allOffset_DVA=[allOffset_DVA, offset_degrees];

allOffset_DVA2=[allOffset_DVA2, ang];

allFlag_rej=[allFlag_rej, tmpRej];
    end
% 
%     [theta, rho] = cart2pol(X, Y);
% polarplot(theta, D, '-o');
% shg
%%

    end% all trials per ppant
    %%
%     clf;


% polarplot(allOffset_theta, allOffset_DVA, 'o');
shg
%% ok, determine range and plot.
%adjust by average height.
% allY = allY- mean(allY);
% allX = allX- mean(allX);
rangeY= [min(allY), max(allY)];
rangeX= [min(allX), max(allX)];

% Yvec = linspace(rangeX(1), rangeX(2), )


% clf
% scatter(allX,allY)

GFX_allTargLocs(ippant).targX = allX;
GFX_allTargLocs(ippant).targY = allY;


GFX_allTargLocs(ippant).eyeX = eyeX;
GFX_allTargLocs(ippant).eyeY = eyeY;


% GFX_allTargLocs(ippant).eyeoffsetDVA = allOffset_DVA;
GFX_allTargLocs(ippant).eyeoffsetDVA = allOffset_DVA2;

GFX_allTargLocs(ippant).eyeoffsetTheta = allOffset_theta;

GFX_allTargLocs(ippant).allFlag_rej = allFlag_rej;



%% append the flag for target removal?
trial_summaryTable.targOutsideDVA = allFlag_rej';


disp(['appending targ-gaze dva offset to trial_summarytable ' filename])
    save(filename, 'trial_summaryTable', '-append');
    

end % ppants
%% 
figure(1);clf
% bug in ppants 18,19.
for ippant= 1:length(GFX_allTargLocs)
subplot(6,7,ippant)
    scatter(GFX_allTargLocs(ippant).targX,GFX_allTargLocs(ippant).targY); hold on
    ylim([-.35 .35])
    xlim([-.25 .25])
% hold on;
%     scatter(GFX_allTargLocs(ippant).targX,GFX_allTargLocs(ippant).targY); hold on
%%

end
shg
figure(2); clf;
clf
% bug in ppants 18,19.
for ippant= 1:length(GFX_allTargLocs)
subplot(6,7,ippant)
    scatter(GFX_allTargLocs(ippant).eyeX,GFX_allTargLocs(ippant).eyeY); hold on
    ylim([-.35 .35])
    xlim([-.25 .25])


end
shg

%%
figure(3); clf;
clf
% bug in ppants 18,19.
allDs_outlier=[];
for ippant= 1:length(GFX_allTargLocs)
subplot(6,7,ippant)
%     scatter(GFX_allTargLocs(ippant).eyeX,GFX_allTargLocs(ippant).eyeY); hold on
 polarplot(GFX_allTargLocs(ippant).eyeoffsetTheta, GFX_allTargLocs(ippant).eyeoffsetDVA, 'o');
% rlim([0 10])


allDs_outlier(ippant,:) =  sum(GFX_allTargLocs(ippant).allFlag_rej);
end
shg

% find()

%% MS figure:
figure(2); clf;

mX=[];

for ippant= 1:length(GFX_allTargLocs)%[1:17,20:34,36:38]

    subplot(2,1,1)
    sc=scatter(GFX_allTargLocs(ippant).targX,GFX_allTargLocs(ippant).targY, 1, 'b'); hold on
xlim([-.18 .18])  
ylim([-.11 .11]);
%   axis equal;
%   hold on;
title('Target pos at onset')


    subplot(212)
        sc=scatter(GFX_allTargLocs(ippant).eyeX,GFX_allTargLocs(ippant).eyeY, 1, 'b'); hold on
xlim([-.18 .18])  
ylim([-.11 .11]);
%   axis equal;
%   hold on;
% mX(ippant)= mean(GFX_allTargLocs(ippant).targX);
end
%
% subplot(211)
%   rectpos = [x,y, w,ht];
rectpos= [-.175, -.1, .35, .2];
rectangle('pos', rectpos);
subplot(212)
rectangle('pos', rectpos);

xlim([-.25 .25])  
ylim([-.15 .15]);
title('Gaze at target onset')
shg
%%
% plot density version
XGrid = [-.175:.0025:.175]; % 35 cm width
YGrid = [-.1:.0025:.1]; % 20 cm height.
densGrid= zeros(length(XGrid), length(YGrid));

% stack each occurence.
for ippant= 1:length(GFX_allTargLocs)%[1:17,20:34,36:38]
   % find nearest X and Y
   assignX = dsearchn(XGrid', GFX_allTargLocs(ippant).targX);
   assignY = dsearchn(YGrid', GFX_allTargLocs(ippant).targY);
   for isamp = 1:length(assignX)
                %increment by 1
                densGrid(assignX(isamp), assignY(isamp)) =   densGrid(assignX(isamp), assignY(isamp)) + 1;
   end

end
%
% densGrid(densGrid==0)=nan;
imagesc(XGrid,YGrid,densGrid);
c=colorbar;
ylabel(c, 'Frequency');
axis equal;
ylabel('Vertical position  [m]');
xlabel('Horizontal position [m]');
set(gca,'fontsize',15)
cmap = colormap;
cmap(1,:) = [.9,.9,.9];
colormap(cmap);
rh=rectangle('pos', rectpos);
rh.LineWidth= 1;
xlim([-.25 .25])  
ylim([-.15 .15]);
set(gcf,'color', 'w'); %axis off