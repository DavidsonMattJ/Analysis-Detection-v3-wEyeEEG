% quick SSVEP check.

%% %events 49 and 50 appear to be start and end.

EEGepochs=[];
epochcount=1;
for iEvent = 1:length(EEG.event)

    if strcmp(EEG.event(iEvent).type, '49') % onset.

        % epoch data.
        onset = EEG.event(iEvent).latency;
        offset = EEG.event(iEvent+1).latency;

        EEGepochs(:,:,epochcount) = EEG.data(:, onset:onset+3000); % 10 sec
    epochcount=epochcount+1;
    end
end

% channels are 14 and 15 for O1 and O2.

% per trial, take the mean pwelch per channels and save.

% spectrum

%% 
usechan=[14,15]
specout= [];
for iepoch = 1:size(EEGepochs,3)

% take the pwelch per trial.
for ichan= 1:2

    tmp = squeeze(EEGepochs(usechan(ichan),:,iepoch));
    [p,freqz]= pwelch(tmp, [],[],[],300);
    specout(ichan,iepoch,:)= log10(p);

    %snr - divide

end


end
% plot mean of epochs.
tmp = mean(squeeze(mean(specout,1)));
subplot(211);
plot(freqz, tmp)
xlim([0 30])

snr1 = conv(tmp, [-.25 -.25  0 0 1 0 0 -.25 -.25], 'same');

hold on;
subplot(212)
plot(freqz,snr1, 'r-')
xlim([0 30])
shg


%% whole recording?
[p,f] = pwelch(EEG.data(14,:),[],[],[],300);
[p2,f] = pwelch(EEG.data(15,:),[],[],[],300);

p3 =(p+p2)/2
snr1 = conv(log10(p), [-.10 -.10 -.10  -.10  -.10  0 0 0 1 0 0 0 -.10 -.10 -.10 -.10  -.10 ], 'same');
snr2 = conv(log10(p2), [-.10 -.10 -.10  -.10  -.10  0 0 0 1 0 0 0 -.10 -.10 -.10 -.10  -.10 ], 'same');
snr3= conv(log10(p3), [-.10 -.10 -.10  -.10  -.10  0 0 0 1 0 0 0 -.10 -.10 -.10 -.10  -.10 ], 'same');

m= (snr1 + snr2)/2
clf
plot(f,  log10(p3) )
shg
%%
tapers=1;
timewin= max(EEG.times)/1000;
freq= f;
noisebinwidth= 2; % Hz to skip
checkshape=1;
etrial = log10(p3)
[kernel, Hbwidth] = buildconvSNRkernel(tapers, timewin, freq', noisebinwidth, checkshape, etrial)
%%
hold on;
plot(f, snr3)
shg;
%%
clf

snr = conv(m, [-.10 -.10 -.10  -.10  -.10  0 0 0 1 0 0 0 -.10 -.10 -.10 -.10  -.10 ], 'same');
hold on;
plot(f, snr, 'ro-')

