%small world analysis.

% boot through GFX params and retain best quadratic fit combination.
GFX_params_qntl=zeros(nsubs,5,2); % slope and thresh
GFX_Acc_perQntl= zeros(nsubs, 5, 7); % qntls and contrasts.

dataIN= GFX_signifit_byGait;


fitstoUse= 2;
extractFits= {'', '_sharedX'}; % so if cfg is set to 1, we extract the shared fit for plots).

for ippant=1:nsubs

    for iqntl= 1:5

        extractfield = ['q' num2str(iqntl) extractFits{fitstoUse}];

       
                NumPerQ= dataIN(ippant, iqntl).(['data_' extractfield])(:,2);
                TotalPerQ= dataIN(ippant, iqntl).(['data_' extractfield])(:,3);
                StimList= dataIN(ippant, iqntl).(['data_' extractfield])(:,1);

        mAcc= NumPerQ./TotalPerQ;

        % thresh, slope.
                dt=dataIN(ippant,iqntl).(['fitresult_' extractfield]);
        %                 dt.Fit
     
        GFX_params_qntl(ippant, iqntl,:)= [dt(1), dt(6)]; % thresh and slope (thresh at 50%).
        GFX_Acc_perQntl(ippant, iqntl,:)= mAcc;
    end % per qntl
end % ippant.


%% 
testmax= fliplr(21:22);%0:24)%30);
bestFit=[];
pvalue=1;
betaM= 0; % we want the maximum negative beta.
tic
for imax= 1:length(testmax)
            thismax= testmax(imax);

    disp(['testing max ' num2str(thismax)]) 
    for iperm= 1:1000
        %select random from list
        thisPerm = randperm(36, thismax);

        % test fit.
         fout = testGoF(GFX_params_qntl(thisPerm,:,2));
        
         % if the correct shape.
         betaV= fout.MDLs.quadMDL.fixedEffects;    
         %we want 
         if fout.compBase_v_Quad.pValue(2)<.05 && fout.compLinear_v_Quad.pValue(2) < .05 &&  betaV(3)<betaM  % negative quadratic
             bestFit.plist = thisPerm;
             pvalue = fout.compBase_v_Quad.pValue(2);
             disp(['Update! ' num2str(pvalue) ',beta = ' num2str(betaV(3))]);
             betaM = betaV(3);
         end
    end
    toc
end
disp('FInished');
toc