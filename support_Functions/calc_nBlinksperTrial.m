% calc_nBlinksperTrial
% - run after J1a_cleanEyeMovementpertrial

setmydirs_detectv3;

cd(procdatadir)
% show ppant numbers:
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%
Fs=90;

interpwindow = 6; % iw * 11ms  @ (90 Hz) % this is in smaples.
% interpwindow = .20; % Sec before and after blink to define interpolation.
%%
GFX_blinksAvpertrial=[];
GFX_blinksoverlappertrial=[];
for ippant = 1:length(pfols)
    cd(procdatadir)


    %     pkdist = participantstepwidths(ippant);
    %%load data from import job.
    load(pfols(ippant).name, 'HeadPos', 'EyePos', 'EyeDir','trialInfo', 'trial_summaryTable', 'subjID');
    savename = pfols(ippant).name;
    disp(['Preparing blink count ' subjID]);

PFX_blinkspertrial=[];
PFX_blinksoverlappertrial=[];
for itrial= 1:length(EyePos)


PFX_blinkspertrial(itrial) = length(EyePos(itrial).blinksAt);

blinkOns = EyePos(itrial).blinksAt;
targOns =find(diff(trialInfo(itrial).targstate)==1);
% targets within 400 ms?
% 400 ms at 90 Hz = 36 samps
 tcount=0; %
if ~isempty(blinkOns)
   
    for itarg = 1:length(targOns) % calc min dist to closest blink.
        minITI = blinkOns - targOns(itarg);

        if min(abs(minITI)<=36) % count if within threshold
            tcount=tcount+1;
        end
    end
    
       
end
PFX_blinksoverlappertrial(itrial) = tcount;
end % trials
GFX_blinksAvpertrial(ippant) = mean(PFX_blinkspertrial);
GFX_blinksoverlappertrial(ippant) = mean(PFX_blinksoverlappertrial);

end

%%
disp(['M blinks per trial = ' sprintf('%.3f',mean(GFX_blinksAvpertrial))])
disp(['SD blinks per trial = ' sprintf('%.3f',std(GFX_blinksAvpertrial))])
disp(['range  = ' sprintf('%.3f',min(GFX_blinksAvpertrial)) '-' sprintf('%.3f',max(GFX_blinksAvpertrial))])


%%
disp(['-------------']);
disp(['M blinks overlap per trial = ' sprintf('%.3f',mean(GFX_blinksoverlappertrial))])
disp(['SD blinks overlap  per trial = ' sprintf('%.3f',std(GFX_blinksoverlappertrial))])
disp(['range  = ' sprintf('%.3f',min(GFX_blinksoverlappertrial)) '-' sprintf('%.3f',max(GFX_blinksoverlappertrial))])