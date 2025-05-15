% check if there is regular spacing between saccade onsets:


setmydirs_detectv3;

cd(procdatadir);

subjFols= dir([pwd filesep '*summary_data.mat']);

visualiseResults=0;

job=[];
job.calcperPpant=1;
job.plotGFX=1;
%%
if job.calcperPpant==1

GFX_saccSpacing=[];

for ippant=1:length(subjFols)
  cd(procdatadir)
  load(subjFols(ippant).name, 'trial_summaryTable', 'EyeDir');

  statTrials = find(trial_summaryTable.isStationary==1);
  statTrIndx = unique(trial_summaryTable.trial(statTrials));

  wlkTrials = find(trial_summaryTable.isStationary==0);
  wlkTrIndx = unique(trial_summaryTable.trial(wlkTrials));
  subjID= subjFols(ippant).name;
  usetrials = {statTrIndx, wlkTrIndx};
    for icond=1:2 % stationary and static.
        
        % find all the trials either stationary or walking. 
        
        

        ppant_saccDiff=[];
        tmpTrials = usetrials{icond};
        for icondtrial= 1:length(tmpTrials)
            itrial = tmpTrials(icondtrial);
        
            % skip if identified as a poor trial:
            skip=0;
            rejTrials_detectv3
            if skip ==1
                continue;
            end


            %% all sacc onsets this trial:
            all_trialsacc = diff(unique([EyeDir(itrial).saccade_startsX,EyeDir(itrial).saccade_startsY]));

        ppant_saccDiff= [ppant_saccDiff, all_trialsacc];
        end

        
            GFX_saccSpacing(ippant,icond,:)= histcounts(ppant_saccDiff,[1:1000]);
        
    end % both conditions
    disp(['fin sacc hist for ippant ' num2str(ippant)]);
end % each ppant.

end % job calc



if job.plotGFX==1
%% plot:
speedCols={'b',[1, 171/255, 64/255], 'k'}; % blue, "normal" yellow

% each sample was 11 ms in unity. (90 Hz).
xvec = [1:size(GFX_saccSpacing,3)]./90;
clf;

box on
for idata=1:2
    tmpD= squeeze(GFX_saccSpacing(:,idata,:));
   
maxperSubj= max(tmpD,[],2);
normGFX= tmpD ./ (repmat(maxperSubj, [1, size(tmpD,2)]));

normGFX = tmpD;
hold on; 
subplot(1,2,idata)

plot(xvec, nanmean(normGFX,1), 'color', speedCols{idata}, 'LineWidth',2); hold on;
bar(xvec, nanmean(normGFX,1), 'FaceColor', speedCols{idata}, 'FaceAlpha', .5); hold on;
errorbar(xvec, nanmean(normGFX,1), CousineauSEM(normGFX),'LineStyle','none', 'Color',[.2 .2 .2])

% bar(xvec, nanmean(normGFX,1), 'FaceColor', speedCols{idata}, 'FaceAlpha', .5)
% errorbar(xvec, nanmean(normGFX,1), CousineauSEM(normGFX),'LineStyle','none', 'Color',[.2 .2 .2])
xlim([0 1]);
% ylim([0 1]);
xlabel('Time since last saccade')
ylabel(' Saccade count')
% 

end
shg
end
%%
% %%
% data =squeeze(mean(normGFX,1)); % your time-series data here %;
% sample_rate =90; % your sample rate here %;
% 
% % Compute the FFT
% N = length(data); % Length of the time-series
% T = 1/sample_rate; % Sampling period
% frequencies = (0:N-1) / (N*T); % Frequency axis
% fft_data = fft(data);
% 
% % Plot the frequency spectrum
% figure(2); clf;
% %plot pos half:
%  plot(frequencies(1:N/2), log(abs(fft_data(1:N/2)))) 
% % plot(frequencies, abs(fft_data));
% xlabel('Frequency (Hz)');
% ylabel('Magnitude');
% title('Frequency Spectrum');



%%%%% FUNCTIONS CALLED
function [fitted,resnorm,t,yInterp] = myDHOfit(y)

% interpd ver:
interpAt = linspace(1, length(y), 1000);
yInterp= interp1(1:length(y), y, interpAt);

yInterp= detrend(yInterp);

% plot(interpAt, propInterp, 'o-', 'Color', pc.Color);

t = interpAt;
% y = exp(-0.5*t).*sin(2*pi*5*t);
oscillator=@(x,t)x(1)*exp(-x(2)*t).*sin(x(3)*t +x(4));

% Define initial parameter values
% amp, damp coeffic, freq, phase shift 
% x0 = [1, .1, 10, 90];

x0 = [.9, .01, 2 ,90];

% Fit the oscillator equation to the data
[fitted,resnorm] = lsqcurvefit(oscillator, x0, interpAt, yInterp);

% return fitted, resnorm, t;

end

function yhat= dampedHarmonic_oscillator(x0,t)
%  amp, damp coeffic, freq, phase shift 
% Define oscillator equation
yhat= x0(1)*exp(-x0(2)*t).*sin(x0(3)*t +x0(4));

end
