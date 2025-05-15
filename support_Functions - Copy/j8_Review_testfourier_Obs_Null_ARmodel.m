
% j8_Review_testfourier_Obs_Null_ARmodel
%
%This job calculates the fourier fits for our data, based on the surrogate
%datasets simulated with an autoregressive model.


setmydirs_detectv3;

%% show ppant numbers:
cd(procdatadir)
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%

%% now with a new shuffling procedure as per Reviewer 3.
nPerm=1000;

%% There are various data types we will want to test
cd([procdatadir filesep 'GFX'])
load('GFX_Data_inGaits.mat');
load('GFX_Data_inGaits_FourierFits.mat') ; % we will append to the presaved, 
%in case completing at separate times.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   There are various datatypes, gaits, etc to step through (and save).
% GFX_FourierNull=[];

for testtype= 3%[1,3,4,5]

    usebin=1; % USE BINNED VERSIONS OF DATA (for now).
%     usebin=0; % use raw to simulate AR model.
    
    switch testtype
        case 1
            %test RT relative to target onset:
            dataIN = GFX_TargPosData;
            typeOnset = 'Target';
            typeDV = 'RT';
        case 2
            %%
            %test RT relative to response onset:
            dataIN = GFX_RespPosData;
            typeOnset = 'Response';
            typeDV = 'RT';
            %%
        case 3
            %test ACC relative to target onset.
            
            dataIN = GFX_TargPosData;
            typeOnset = 'Target';
            typeDV = 'Accuracy';
            
        case 4
            % Test response (click) likelihood, relative to resp onset.
            dataIN = GFX_RespPosData;
            typeOnset = 'Response';
            typeDV='Counts';
        case 5
            % Test response (click) likelihood, relative to resp onset.
            dataIN = GFX_TargPosData;
            typeOnset = 'Target';
            typeDV='Counts';
            
        case 6
             dataIN = GFX_SaccadePosData;
            typeOnset = 'Saccade';
            typeDV='Counts';
    end
    %%
    cfg=[];
    cfg.subjIDs = subjIDs;
    cfg.type = typeOnset;
    cfg.DV = typeDV;
%     cfg.datadir= datadir; % for orienting to figures folder
    cfg.HeadData= GFX_headY;
    cfg.pidx1= pidx1;
    cfg.pidx2= pidx2;
    cfg.plotlevel = 'GFX'; % plot separate figures per participant
    cfg.norm=0; % already z scored, so don't tweak.
    cfg.ylims = [-.15 .15]; % if norm =0;
    cfg.normtype= 'relative';
    %%
    
    %just one gait at a time (its a slow process).
   for  nGaits_toPlot=2%1:2
    
    % plot_FourierFit(cfg,dataIN);
    %%
   
   
        iLR=3;
    gaitfield = {'gc', 'doubgc'};
    binfield = {'','_binned'};
    
    ppantData=[];
    
    shuffData=[];
    
    %which field of datastructure to plot?
    if strcmp(cfg.DV, 'RT')
        usefield = [gaitfield{nGaits_toPlot} binfield{usebin+1} '_rts'];
        ylabis = 'z(RT)';
    elseif strcmp(cfg.DV, 'Accuracy')
        usefield = [gaitfield{nGaits_toPlot} '_binned_Acc'];
        if ~cfg.norm
            ylabis=  cfg.DV;
        else
            ylabis = [cfg.DV 'norm: ' cfg.normtype];
        end
    elseif strcmp(cfg.DV, 'Counts')
        usefield = [gaitfield{nGaits_toPlot} '_binned_counts'];
        ylabis = [cfg.type ' ' cfg.DV];

        if strcmp(cfg.type, 'Saccade')
 usefield = [gaitfield{nGaits_toPlot} '_sacc_all_binned_counts'];

        end
    end
    
    %collate data:
    allAR_models=[];
    nppants= size(dataIN,1);
    for isub= 1:size(dataIN,1)
        
        ppantData(isub,:)= dataIN(isub,iLR).(usefield);
        
        %% while here collect the model fits:
        
ar_model = ar(ppantData(isub,:), 1);

% Get the AR coefficients
ar_coeffs = ar_model.A;
        
     allAR_models(isub,:) = ar_coeffs;
    %%
    end
   
    
    %% extract fourier fits per shuffled series:
    % Declaring the type of fit.
    FitType = 'fourier1';
    % Creating and showing a table array to specify bounds.
    CoeffNames = coeffnames(fittype(FitType));
    %%
    %set bounds for w
    CoeffBounds = array2table([-Inf(1,length(CoeffNames));...
        Inf(1,length(CoeffNames))],'RowNames',...
        ["lower bound", "upper bound"],'VariableNames',CoeffNames);
    %%
    % Specifying bounds according to the position shown by the table.
    % e.g. to force fit with w ~ 1.545, we ChangeBound of w parameter
    % CoeffBounds.w(1) = 1.54;
    % CoeffBounds.w(2) = 1.55;
    %
    Hzspace = [0.01:.2:10];
    % perW=per
    fits_Rsquared_obsrvd = nan(1, length(Hzspace));
%     if fitShuffData==1
        fits_Rsquared_shuff = nan(nPerm, length(Hzspace));
        fits_Rsquared_shuffCV = nan(3, length(Hzspace));
%         meanShuff = squeeze(mean(shuffData,1));
%     end
    
    % for a bunch of perms,     
    % 1) create surrogate AR data. (per subject)
    % 2)take binned avg as per (j4__)
    % 3) take mean over subjs (each perm).
    % 4) fit all freqs, each perm GFX.
    
      pidx=ceil(linspace(1,100,41));% 
      
    groupData = squeeze(mean(ppantData));
      
    
    gM= groupData;
    
    
    %x axis:          %approx centre point of the binns.
    mdiff = round(mean(diff(pidx)./2));
    xvec = pidx(1:end-1) + mdiff;
    %%
    

        %%

    disp(['creating all surrogates...']);
    storetestData=zeros(nPerm,length(gM));

    
    nPerm=1000;
%     it an AR model with order 1
order = 1;
ar_model = ar(data, 1);

% Get the AR coefficients
ar_coeffs = ar_model.A;
n=length(gM);

mean_empirical = mean(gM);
std_empirical = std(gM);

% lets make fake group data, based on individual AR fits.
% this way (realising new group data based on sub-level AR surrogates,
%is preferable according to Re et al., 2022).

%% first create a bank of fake ppant AR data.
% this method create data per ppant:
% nfake= 50;
% allppant_ARbank= nan(nppants, nfake ,n);
% 
% for ippant = 1:size(ppantData,1)
%    % retreive model fit, spam 10?
%    ppant_fake=nan(nfake, n);  
%    ppant_fit = allAR_models(ippant,:);
%    mean_empirical = mean(ppantData(ippant,:));
%    std_empirical = std(ppantData(ippant,:));
%    
%    for ifake = 1:nfake
%        
%        %% make a fake trial.
%         % Generate white noise for the surrogate data
%         white_noise = randn(n, 1);
%         
%         % create 
%         
%         % Initialize the surrogate data
%         surrogate_data = zeros(n, 1);
% 
%         for i = order + 1:n
%             for j = 1:order
%                 surrogate_data(i) = surrogate_data(i) + ppant_fit(j) * surrogate_data(i - j);
%             end
%             surrogate_data(i) = surrogate_data(i) + white_noise(i);
%         end
%         
%         surrogate_data= detrend(surrogate_data);
%         
%         % Calculate the mean and standard deviation of the surrogate data
%         mean_surrogate = mean(surrogate_data);
%         std_surrogate = std(surrogate_data);
%         
%         % Adjust the surrogate data to match the mean and standard deviation of the empirical data
%         surrogate_data2 = (surrogate_data - mean_surrogate) * (std_empirical / std_surrogate) + mean_empirical;
%         %  Adjust the surrogate data to match the mean of the empirical data
%         % surrogate_data2 = surrogate_data + (mean_empirical - mean_surrogate);
% %         storetestData(iperm,:) = surrogate_data2';
%        
%        allppant_ARbank(ippant,ifake,:) = surrogate_data2;
%        
%        
%    end
%     
% end
% disp(['done'])

%% this method produces nPerm baed on group data.
    for iperm= 1:nPerm


        
        
        % Generate white noise for the surrogate data
        white_noise = randn(n, 1);
        
        % create 
        
        % Initialize the surrogate data
        surrogate_data = zeros(length(gM), 1);

        for i = order + 1:n
            for j = 1:order
                surrogate_data(i) = surrogate_data(i) + ar_coeffs(j) * surrogate_data(i - j);
            end
            surrogate_data(i) = surrogate_data(i) + white_noise(i);
        end

        surrogate_data= detrend(surrogate_data);
% Calculate the mean and standard deviation of the surrogate data
mean_surrogate = mean(surrogate_data);
std_surrogate = std(surrogate_data);

% Adjust the surrogate data to match the mean and standard deviation of the empirical data
surrogate_data2 = (surrogate_data - mean_surrogate) * (std_empirical / std_surrogate) + mean_empirical;
%  Adjust the surrogate data to match the mean of the empirical data
% surrogate_data2 = surrogate_data + (mean_empirical - mean_surrogate);
storetestData(iperm,:) = surrogate_data2';
    end
disp(['done!']);




% now instead, our new group data is a random trial from each participant:
   % for each shuffle as well, calc the fits.
for iperm=1:nPerm 
    
    % each perm will use a unique surrogate group level data
    %create new group level result
%         trialPerppant = randi(nfake,[1,nppants]);
%     
%         fake_GM=[];
%         for ippant = 1:size(allppant_ARbank,1);
%            fake_GM(ippant,:) = squeeze(allppant_ARbank(ippant,trialPerppant(ippant), :)); 
%             
%         end
% 
%         test_gM = mean(fake_GM);
%         
       %%
       
    % step through w, forcing fit at particular periods, by updating the
    % bounds in fit options.
    for ifreq= 1:length(Hzspace)
        % include period and Rsquared
        %treat max xvec as our full 'period'
        %             Hzapp = xvec(end)/ (2*pi/(f.w));
        testw = 2*pi*Hzspace(ifreq)/xvec(end);
        
        CoeffBounds.w(1) = testw;
        CoeffBounds.w(2) = testw;
        
        %update fit opts settings
        
        %Update Fit Options setting.
        FitOpts = fitoptions('Method','NonlinearLeastSquares','Lower',table2array(CoeffBounds(1,:)),...
            'Upper',table2array(CoeffBounds(2,:)));
        
        %first test this period on observed data
        %set last coefficient value in the fourier model (w) to this value:
        if iperm==1
        [tmpf,gof] = fit(xvec', gM', 'fourier1', FitOpts);
        % how good/bad was the fit?
        fits_Rsquared_obsrvd(1,ifreq) = gof.rsquare;
        end
    
            try
            [tmpf,gof] = fit(xvec', storetestData(iperm,:)', 'fourier1', FitOpts);
            
%             [tmpf,gof] = fit(xvec', test_gM', 'fourier1', FitOpts);
            
            fits_Rsquared_shuff(iperm,ifreq) = gof.rsquare;
            catch
                disp([' skipping perm ' num2str(iperm) ' ' cfg.type ' ' cfg.DV ...
            ', gaits(' num2str(nGaits_toPlot) ') likely a NaN'])
            end

            
            
    end %ifreq

    disp(['fin perm ' num2str(iperm)]);
      
    end % per perm
    
    
    %quick save?
    %% for visualization, show some individual AR surroates and spectra, as well as their averages.
    % 10 random examples.
    %566, 629, 93, 206
%     showt = randi(1000,3);
% clf;
%     for itrial = 1:3
%     subplot(3,2, 1 + 2*(itrial-1))
%     plot(storetestData(showt(itrial),:)); 
%     subplot(3,2, 2 + 2*(itrial-1))
%    plot(Hzspace, fits_Rsquared_shuff(showt(itrial),:)); 
%         title(num2str(showt(itrial)));
%     end
%     shg
    
    %
    
    %% now we take the percentiles based on the max over all freqs (each shuffle).
    
     % per freq, store the 95%CI. of Rsq values(plotted below)
     
     for ifreq=1:size(fits_Rsquared_shuff,2)   
     fits_Rsquared_shuffCV(:,ifreq) = quantile(fits_Rsquared_shuff(1:nPerm,ifreq), [.05, .5, .95]);
     end
     %find max per perm.
     maxperPerm = max(fits_Rsquared_shuff,[],2); % max at all freqs  
     fits_Rsquared_shuffCV_new = quantile(maxperPerm(1:nPerm), [.05, .5, .95]);
    
    %store:
    GFX_FourierNull.([cfg.type 'Ons_' usefield '_fitsRsq_Obs']) = fits_Rsquared_obsrvd;
    
    GFX_FourierNull.([cfg.type 'Ons_' usefield '_fitsRsq_ShuffCV_AR']) = fits_Rsquared_shuffCV;
    GFX_FourierNull.([cfg.type 'Ons_' usefield '_fitsRsq_ShuffCV_max_AR']) = fits_Rsquared_shuffCV_new;
    
    
    %% REVIEW FIGURE! 
    % sanity checks
    
%     
%     clf;
%     subplot(2,3,3);
%     % overlay
%     showt= [566, 629, 93, 206];
%     for itrial= 1:4
%         subplot(2,3,2);
%         offset= storetestData(showt(itrial),:) - mean(storetestData(showt(itrial),:)) + itrial*.05;
%         
%         plot(offset, 'linew',2);
%         xlabel('stride completion')
%         set(gca,'xtick', [], 'ytick', []);
%         ylabel('a.u.')
%         title('example AR trials')
%         set(gca,'fontsize',15);
%         hold on
%         subplot(2,3,5); hold on
%         plot(Hzspace, fits_Rsquared_shuff(showt(itrial),:), 'linew',2);
%         xlabel('Frequency (cps)');
%         title('(fits)')
%     end
%     set(gca,'fontsize',15);
% %     clf
%     subplot(131)
% figure(itype);
%     oh=plot(Hzspace, fits_Rsquared_obsrvd, 'b', 'linew',2);
%     hold on
% %     plot(Hzspace, fits_Rsquared_shuffCV(2,:), '-k','linew',2);
% %     plot(Hzspace, fits_Rsquared_shuffCV(1,:), ':k','linew',2);
%     u95h=plot(Hzspace, fits_Rsquared_shuffCV(3,:), ':k','linew',2);
%     u95h2=plot(xlim, [fits_Rsquared_shuffCV_new(3),fits_Rsquared_shuffCV_new(3)], ':r', 'linew',2);
%     hold on;
% %     plot(Hzspace, fits_Rsquared_shuffAv, ':m','linew',2);
%     shg
%     shg
%     legend([oh, u95h, u95h2],...
%         {'Accuracy (empirical)', '95% CI (aperiodic null)', '95% CI (all freq.)'});
%     title('Aperiodic surrogate')
%     xlabel('Frequency (cps)');
%     set(gca,'fontsize',15);
%     ylim([0 1]);
%     %% show some example trials
%     
%     %%
%     % hold and plot OG for comparison.
%     %% plot the average result:
%     subplot(133)
%       oh=plot(Hzspace, fits_Rsquared_obsrvd, 'b', 'linew',2);
%     hold on
%     % for 1000 reps, take average of N=35 surrogates.
%     pAR = [];
%     fits_Rsquared_shuffCVAR_mean=[];
%     for ip= 1:nPerm;
%         
%         avover = randi(1000, [1,35]);
%         pAR(ip,:) = mean(fits_Rsquared_shuff(avover,:),1);
%     end
%     % plot meana nd quantiles:
%     for ifreq= 1:size(pAR,2)
%     fits_Rsquared_shuffCVAR_mean(:,ifreq) = quantile(pAR(:,ifreq), [.05, .5, .95]);
%     end
%     title('Alternate aperiodic surrogate') 
%     u95h=plot(Hzspace, fits_Rsquared_shuffCVAR_mean(3,:), ':k','linew',2);
%     %find max per perm.
%      maxperPerm = max(pAR,[],2); % max at all freqs  
%      fits_Rsquared_shuffCVAR_new = quantile(maxperPerm(1:nPerm), [.05, .5, .95]);
%      u95h2=plot(xlim, [fits_Rsquared_shuffCVAR_new(3),fits_Rsquared_shuffCVAR_new(3)], ':r', 'linew',2);
%     hold on;
%      legend([oh, u95h, u95h2],...
%         {'Accuracy (empirical)', '95% CI (mean spectra)', '95% CI (all freq.)'});
%     xlabel('Frequency (cps)');
%     set(gca,'fontsize',15);
%     ylim([0 1]);
    %%
    
    disp(['finished gait ' num2str(nGaits_toPlot) ', ' typeOnset ' ' typeDV])
end % nGaits

end %itype
GFX_FourierNull.nShuff = size(meanShuff,1);
%%
% save('GFX_Data_inGaits_FourierFits', 'GFX_FourierNull', 'Hzspace')
save('GFX_Data_inGaits_FourierFits', 'GFX_FourierNull', 'Hzspace','-append')