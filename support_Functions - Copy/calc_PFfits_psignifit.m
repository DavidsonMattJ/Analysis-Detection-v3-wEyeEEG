function calc_PFfits_psignifit(dataIN, cfg)


% as an alternative to the palamedes toolbox, here trying the psignifit
% function library. More flexibility and parameters in the results.

% dataIN= GFX_TargPosData;
%pseudo:
%1) load ind data in the format of a matrix (x | nCorr | nTotal) .
%2) Fit per ppant and plot various types of fits.

usecolsStWlk = {[.7 .7 .7], [.7 .7 0]}; % R Gr
usecolsLR= {[.7 0, 0], [0 .7 0], [.7, 0, .7]}; % R Gr Prple.
figure(1); clf; set(gcf, 'color', 'w', 'units', 'normalized', 'position', [0 0 .9 .9]);
nGaits_toPlot=1; % no double dipping.

usegaitfields = {'gc_', 'doubgc_'};
use_conds= {'_allstatnry', '_allwlking'};
use_steps = {'_LRwlking'};
           
useg= usegaitfields{nGaits_toPlot};
 
useLogscale=0;


%% fit options:
options =[];
options.sigmoidName = 'norm';   % choose a cumulative Gaussian as the sigmoid
options.expType= 'YesNo'; % free params for lapse and guess rate.
%%
% % another standard alternative is the logistic function
% options.sigmoidName    = 'logistic';
% 
% % For data on a logscale you may want to fit a log-normal distribution or a
% % Weibull which you invoce with:
% options.sigmoidName    = 'logn';
% %or
% options.sigmoidName    = 'weibull';
% 
% % We also included the gumbel and reversed gumbel functions for asymmetric
% % psychometric functions. The gumbel has a longer lower tail the reversed
% % gumbel a longer upper tail. 
% 
% options.sigmoidName    = 'gumbel';
% % or
% options.sigmoidName    = 'rgumbel';
% 
% % for a heavy tailed distribution use
% options.sigmoidName    = 'tdist';
% 

% %preallocate for storage:
contrVals=[];
PFX_Fits_StWlk_xy=[];

job.calcconcat_pfx=1;

if job.calcconcat_pfx==1

    % the concat job pre fits for 
    % all standing and walking
    % walking (left or right foot).
    % all walking (byquintike).
    GFX_signifit=[];
    GFX_signifit_byGait=[];

for ippant=1:size(dataIN,1)
    
    % first extrac thte data for walking, single steps, different feet:
  disp(['Performing psignif fits for ppant ' num2str(ippant) '/ ' num2str(size(dataIN,1))]);
            %%
            for iLR=1:2 % left / right step leading.
                

                NumPer = dataIN(ippant,iLR).([useg 'NumPerContr_LRwlking']);
                TotalPer = dataIN(ippant,iLR).([useg 'TotalPerContr_LRwlking']);
                StimList= dataIN(ippant,iLR).(['gc_ContrVals_LRwlking']);
                % add to output:

                RTs = dataIN(ippant,iLR).(['gc_RTperContr_LRwlking']);


                % adjust stimList to show diff between background and
                % target?
%                 StimList= StimList-0.4; % .4 was bg Grey.

                dataFit = [StimList', NumPer', TotalPer'];

                % perform fit:
                fitresult = psignifit(dataFit,options);
                resultSmall = rmfield(fitresult,{'Posterior','weight'});
                clear fitresult % keep mem free as we go
                % fake plot!
                set(gcf,'visible', 'off'); % hide figure, but we want the output:
                pr=plotPsych(resultSmall);
                %result.fit = T W L G;
                %     slope = (1 - G - L) / (4 * W) **( at threshold)

                slope = (1 - resultSmall.Fit(4) - resultSmall.Fit(3)) / (4* resultSmall.Fit(2));
    
                GFX_signifit(ippant,iLR).(['lowX_LRwlking'])= pr.low.XData;
                GFX_signifit(ippant,iLR).(['lowY_LRwlking' ])= pr.low.YData;
                GFX_signifit(ippant,iLR).(['mainX_LRwlking'])= pr.main.XData;
                GFX_signifit(ippant,iLR).(['mainY_LRwlking'])= pr.main.YData;
                GFX_signifit(ippant,iLR).(['highX_LRwlking'])= pr.high.XData;
                GFX_signifit(ippant,iLR).(['highY_LRwlking'])= pr.high.YData;
                GFX_signifit(ippant,iLR).(['fitresult_LRwlking'])= [resultSmall.Fit', slope];
                GFX_signifit(ippant,iLR).(['data_LRwlking'])= [dataFit,RTs'];


            end % i LR
         
% now combined case, both feet together:
iLR=3;
for iStWlk=1:2

    %note third index is all single steps combined:
    NumPer = dataIN(ippant,3).([useg 'NumPerContr' use_conds{iStWlk}]);
    TotalPer = dataIN(ippant,3).([useg 'TotalPerContr' use_conds{iStWlk}]);
    StimList= dataIN(ippant,3).(['gc_ContrVals' use_conds{iStWlk}]);

    
                % adjust stimList to show diff between background and
                % target?
%                 StimList= StimList-0.4; % .4 was bg Grey.


     RTs = dataIN(ippant,3).(['gc_RTperContr' use_conds{iStWlk}]);
    % wrangle data:
    dataFit = [StimList', NumPer', TotalPer'];


    %%
    fitresult = psignifit(dataFit,options);
    resultSmall = rmfield(fitresult,{'Posterior','weight'});
     clear fitresult % keep mem free as we go
    % plot!    
    set(gcf,'visible', 'off'); % hide figure, but we want the output:
    pr=plotPsych(resultSmall);

    %result.fit = T W L G;    
%     slope = (1 - G - L) / (4 * W) **( at threshold)
    slope = (1 - resultSmall.Fit(4) - resultSmall.Fit(3)) / (4* resultSmall.Fit(2));
    
    GFX_signifit(ippant,3).(['lowX' use_conds{iStWlk}])= pr.low.XData;
    GFX_signifit(ippant,3).(['lowY' use_conds{iStWlk}])= pr.low.YData;
    GFX_signifit(ippant,3).(['mainX' use_conds{iStWlk}])= pr.main.XData;
    GFX_signifit(ippant,3).(['mainY' use_conds{iStWlk}])= pr.main.YData;
    GFX_signifit(ippant,3).(['highX' use_conds{iStWlk}])= pr.high.XData;
    GFX_signifit(ippant,3).(['highY' use_conds{iStWlk}])= pr.high.YData;
    GFX_signifit(ippant,3).(['fitresult' use_conds{iStWlk}])= [resultSmall.Fit', slope];
    GFX_signifit(ippant,3).(['data' use_conds{iStWlk}])= [dataFit, RTs'];
    
end % standing or walking.



% for the final case, we will perform 5 fits per gait (quintiles).
NumPerQ= dataIN(ippant, 3).(['gc_NumPerContr_qntlwlking']);
TotalPerQ = dataIN(ippant,3).(['gc_TotalPerContr_qntlwlking']);
StimListQ= dataIN(ippant,3).(['gc_ContrVals_qntlwlking']);
RTsperQ= dataIN(ippant,3).(['gc_RTperContr_qntlwlking']);

for iqnt= 1:size(NumPerQ,1)
    %qntl specific data:
    NumPer = NumPerQ(iqnt,:);
    TotalPer = TotalPerQ(iqnt,:);
    Stimlist = StimListQ(iqnt,:);


                % adjust stimList to show diff between background and
                % target?
%                 Stimlist= Stimlist-0.4; % .4 was bg Grey.

    RTs = RTsperQ(iqnt,:);
    dataFit = [Stimlist', NumPer', TotalPer'];


    %%
    fitresult = psignifit(dataFit,options);
    resultSmall = rmfield(fitresult,{'Posterior','weight'});
    clear fitresult % keep mem free as we go
    % plot!
    set(gcf,'visible', 'off'); % hide figure, but we want the output:
    pr=plotPsych(resultSmall);
    %result.fit = T W L G;    
%     slope = (1 - G - L) / (4 * W) **( at threshold)
    slope = (1 - resultSmall.Fit(4) - resultSmall.Fit(3)) / (4* resultSmall.Fit(2));

    GFX_signifit_byGait(ippant,iqnt).(['lowX_q' num2str(iqnt)])= pr.low.XData;
    GFX_signifit_byGait(ippant,iqnt).(['lowY_q' num2str(iqnt)])= pr.low.YData;
    GFX_signifit_byGait(ippant,iqnt).(['mainX_q' num2str(iqnt)])= pr.main.XData;
    GFX_signifit_byGait(ippant,iqnt).(['mainY_q' num2str(iqnt)])= pr.main.YData;
    GFX_signifit_byGait(ippant,iqnt).(['highX_q' num2str(iqnt)])= pr.high.XData;
    GFX_signifit_byGait(ippant,iqnt).(['highY_q' num2str(iqnt)])= pr.high.YData;
    GFX_signifit_byGait(ippant,iqnt).(['fitresult_q' num2str(iqnt)])= [resultSmall.Fit', slope];
    GFX_signifit_byGait(ippant,iqnt).(['data_q' num2str(iqnt)])= [dataFit, RTs'];


    % also save for easier comparison, the fit with a shared stimlist
    % (average).
    %adjust for BG  
    tmpStim = StimListQ-0.4;

    dataFit = [mean(StimListQ,1)', NumPer', TotalPer'];
 
    fitresult = psignifit(dataFit,options);
    resultSmall = rmfield(fitresult,{'Posterior','weight'});
    clear fitresult % keep mem free as we go
    % plot!
    set(gcf,'visible', 'off'); % hide figure, but we want the output:
    pr=plotPsych(resultSmall);
    %result.fit = T W L G;    
%     slope = (1 - G - L) / (4 * W) **( at threshold)
    slope = (1 - resultSmall.Fit(4) - resultSmall.Fit(3)) / (4* resultSmall.Fit(2));

    GFX_signifit_byGait(ippant,iqnt).(['lowX_q' num2str(iqnt) '_sharedX'])= pr.low.XData;
    GFX_signifit_byGait(ippant,iqnt).(['lowY_q' num2str(iqnt) '_sharedX'])= pr.low.YData;
    GFX_signifit_byGait(ippant,iqnt).(['mainX_q' num2str(iqnt) '_sharedX'])= pr.main.XData;
    GFX_signifit_byGait(ippant,iqnt).(['mainY_q' num2str(iqnt) '_sharedX'])= pr.main.YData;
    GFX_signifit_byGait(ippant,iqnt).(['highX_q' num2str(iqnt) '_sharedX'])= pr.high.XData;
    GFX_signifit_byGait(ippant,iqnt).(['highY_q' num2str(iqnt) '_sharedX'])= pr.high.YData;
    GFX_signifit_byGait(ippant,iqnt).(['fitresult_q' num2str(iqnt) '_sharedX'])= [resultSmall.Fit', slope];
    GFX_signifit_byGait(ippant,iqnt).(['data_q' num2str(iqnt) '_sharedX'])= [dataFit, RTs'];

 slp1=GFX_signifit_byGait(ippant,iqnt).(['fitresult_q' num2str(iqnt)])(6);
 slp2=GFX_signifit_byGait(ippant,iqnt).(['fitresult_q' num2str(iqnt) '_sharedX'])(6);
 
    disp(['qfit slope= ' num2str(slp1) ', vs shared = ' num2str(slp2)]);
end

    end % ppant; 


%% store the GFX structure.
save('GFX_Data_inGaits_SigniFITs.mat', 'GFX_signifit', 'GFX_signifit_byGait'); % compressed format
%%

end % calcconat job:


  % psThresh = fitresult.Fit(1);
    % psWidth = fitresult.Fit(2);
    % result.conf_Intervals

    % This gives you the basic result of your fit. The five values reported are:
    %    the threshold
    %    the width (difference between the 95 and the 5 percent point of the unscaled sigmoid)
    %    lambda, the upper asymptote/lapse rate
    %    gamma, the lower asymptote/guess rate
    %    eta, scaling the extra variance introduced (a value near zero indicates
    %         your data to be basically binomially distributed, whereas values
    %         near one indicate severely overdispersed data)
    % The field conf_Intervals returns credible intervals for the values provided
    % in options.confP. By default these are 68%, 90% and 95%. With default settings
    % you should thus receive a 5x2x3 array, which contains 3 sets of credible intervals
    % (lower and upper end = 2 values) for each of the 5 parameters.


%     NB:
% The psychometric function is a commonly used model in psychophysics to 
% describe the relationship between a stimulus intensity and a subject's 
% response. The slope of the psychometric function is an important 
% parameter that describes the rate of change in the subject's response as
%  a function of the stimulus intensity.
% 
% The psychometric function can be modeled using a cumulative distribution
% function (CDF) such as the logistic function:
% 
% P(x) = G + (1 - G - L) * (1 / (1 + exp(-(x - T) / W)))
% 
% where:
% 
% x is the stimulus intensity
% P(x) is the proportion of trials in which the subject responds "yes"
% G is the guess rate (the proportion of trials in which the subject guesses)
% L is the lapse rate (the proportion of trials in which the subject lapses)
% T is the threshold (the stimulus intensity at which the subject responds
%  "yes" 50% of the time)
% W is the width or slope parameter that determines how steeply the CDF 
% rises from the threshold
% 
% To compute the slope of the psychometric function, we need to take the
%  derivative of the CDF with respect to the stimulus intensity x:
% 
% dP(x) / dx = (1 - G - L) * (exp(-(x - T) / W) / (W * (1 + exp(-(x - T) / W))^2))
% 
% At the threshold T, the derivative of the CDF is maximal and equal to:
% 
% dP(T) / dx = (1 - G - L) / (4 * W)
% 
% Thus, the slope of the psychometric function at the threshold is:
% 
% slope = dP(T) / dx = (1 - G - L) / (4 * W)
% 
% Therefore, to compute the slope of the psychometric function given the threshold T, width W, guess rate G, and lapse rate L, we can use the formula:
% 
% slope = (1 - G - L) / (4 * W)
% 
% Note that the slope is inversely proportional to the width parameter W, so a narrower psychometric function (i.e., a steeper slope) corresponds to a larger value of W.
