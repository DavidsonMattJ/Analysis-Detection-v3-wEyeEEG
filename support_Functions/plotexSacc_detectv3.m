
% function plotexSacc_detectv3(cfg)
%%

setmydirs_detectv3;
cd(procdatadir);
%
figure(2);
clf
% these specs were used pre:review:
% subjID= 'GC';
% pfile = dir([pwd filesep subjID '*' '.mat']);
% fntsize= 10
% nat_trial= 55;
% now (different trial to show blink detection).

subjID= 'ABB';
pfile = dir([pwd filesep subjID '*' '.mat']);
fntsize= 10
nat_trial= 61;

showBlink=1;
%



% slow_trial=9;
load(pfile.name, 'EyeDir','EyePos','trialInfo');
slowEye = EyeDir;
% load([subjID '_naturalWalking_gazeDirection'], 'EyeDir');
naturalEye= EyeDir;
%
% fntsize=12
subplot(2,3,1:2);

% timevec = linspace(0,6.5,585) -1; % 11 ms samples
timevec = trialInfo(nat_trial).times';    


% if showing blink, plot that first (to be overlayed)

if showBlink==1
% show the blinks interpolated on this trial.
hold on
iw= 6; % interpwindow = 6; % iw * 11ms  @ (90 Hz) % this is in smaples.
blinksAt = EyePos(nat_trial).blinksAt;
blinksEnd = EyePos(nat_trial).blinksEnd;
for iblink = 1:length(blinksAt)

    xvec=[blinksAt(iblink), blinksAt(iblink), blinksEnd(iblink), blinksEnd(iblink)];
xvec =timevec(xvec);
yvec=[-.5 .5 .5 -.5];

ph=patch(xvec,yvec, [.2 .2 .2]);
ph.FaceAlpha= .1;
ph.EdgeColor= 'w';


plot(timevec, EyeDir(nat_trial).Z, 'k'); hold on;
plot(timevec, EyeDir(nat_trial).Y, 'k');


% pg.
end
end

% now plot gaze direction



hold on
pY=plot(timevec, naturalEye(nat_trial).Y_clean, 'color','b', 'LineWidth',1); hold on;
pX=plot(timevec, naturalEye(nat_trial).Z_clean, 'color', 'r', 'LineWidth',1); hold on;

% plot(timevec, naturalHead(nat_trial).Y, 'color', speedCols{2}, 'LineWidth',2); hold on;
shg
% xlim([-1 5.5])
set(gca,'units','normalized', 'fontsize',fntsize,'XTickLabel',[], 'YTick', [-1.2:.2:1.2]);

ylabel({['Gaze direction'];['(m)']}); hold on;
axis tight
ylim([-1.4 1.4])
% axis tight
yls = get(gca,'ylim');
tOnsets = find(trialInfo(nat_trial).targstate ==1);
for itarg = 1:length(tOnsets);
    plot([timevec(tOnsets(itarg)) timevec(tOnsets(itarg))],[yls(1), yls(1)+.3], 'k-');
        plot([timevec(tOnsets(itarg)) timevec(tOnsets(itarg))],[yls(1)+.29], 'dk-');
end
text(0, yls(1)+.2, 'Target Onsets', 'Color','k')
%
% legend
lg=legend([pY,pX], {'Vertical','Horizontal'},'Location', 'NorthWest','Orientation','horizontal');
lg.ItemTokenSize(1)=8;

box on

%%


% set(gca, '')
%
curr_ax = get(gca,'position');
%new pos is (bottom left) x, y, width, height 
offsetY = curr_ax(2)-  curr_ax(4)*.6- .001;
offsetX = curr_ax(1); % none
newpos = [offsetX, offsetY, curr_ax(3) curr_ax(4)*.6];

newax = subplot('position' ,newpos);

% plot the derivative
% copy pasted from the saccade detection script: 
%j2C_timestamp_Saccdes



            % where to perform the threshold detection on?
            timevec = trialInfo(nat_trial).times';
            stopAt = length(timevec);
            eyePos=[];
            eyePos(:,1) = EyeDir(nat_trial).Z_clean(1:stopAt); % note we are clipping the final 500 ms.
            eyePos(:,2) = EyeDir(nat_trial).Y_clean(1:stopAt);
            % Calculate velocity time series
            % first deriv
            velX = diff(eyePos(:,1))./diff(timevec(1:stopAt))'; % Horizontal velocity
            velY = diff(eyePos(:,2))./diff(timevec(1:stopAt))'; % Vertical velocity
            velX = movmean(velX, 5); % Moving average over 5 data points
            velY = movmean(velY, 5); % Moving average over 5 data points

            % Compute velocity thresholds
            velX_med = median(velX); % Median of horizontal velocity
            velY_med = median(velY); % Median of vertical velocity
            velX_std = median(abs(velX - velX_med)); % Median estimator of standard deviation of horizontal velocity
            velY_std = median(abs(velY - velY_med)); % Median estimator of standard deviation of vertical velocity
            velX_thresh = velX_med + 6 * velX_std; % Threshold for horizontal velocity
            velY_thresh = velY_med + 6 * velY_std; % Threshold for vertical velocity

            % plot the thresh and crossings;
plot(timevec(1:length(velX)) ,(abs(velX)), 'r');
hold on;
plot(timevec(1:length(velX)),(abs(velY)), 'b');
pv= plot([timevec(1) timevec(length(velX))], [velX_thresh velX_thresh], 'r:','linew',1);
ph=plot([timevec(1) timevec(length(velY))], [velY_thresh velY_thresh], 'b:', 'LineWidth',1);
axis tight;

ylim([ 0 2*max([velX_thresh, velY_thresh])])

ylabel({['Velocity'];['(m/s)']})

% also plot the detected saccades!
saccX= EyeDir(nat_trial).saccade_startsX;
saccY= EyeDir(nat_trial).saccade_startsY;
for ix= 1:length(saccX)

plot([timevec(saccX(ix)), timevec(saccX(ix))], [velX_thresh velX_thresh], ['ro'], 'linew',2 )

end

for iy= 1:length(saccY)

plot([timevec(saccY(iy)), timevec(saccY(iy))], [velY_thresh velY_thresh], ['bo'], 'linew',2 )

end


lg=legend([pv, ph], {'Vert. threshold', 'Horiz. threshold'}, 'orientation', 'horizontal');
lg.ItemTokenSize(1)=8;
xlabel('Trial time (sec)')
% ylim([0 8])
set(gca,'units','normalized', 'fontsize',fntsize);
%% % now plot the saccade location data:
% from j5A_storeEyeDir_Egocentric

