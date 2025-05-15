% j12_splitHalfReliabilityEstimate

% this script performs per participant, the basic fourier fit analysis.
% i.e. for accuracy, rt, and response onsets over the full gaitcycle.

%Importantly, now we check the split-half reliability. 
% perform many splits, and then compute the split half prob of an R^2 value
% at a particular peak.

%so we can assess how likely a peak is at given cycles per degree, per
%participant.


% r_split_half = 2 * r / (1 + r);


%psuedo-code:
% load participant level table
% 1)subselect half the trials 
% 2)- perform doubgc fit at 2 and 4 Hz. retain R2
% 3)-store per split
% 4) repeat 2-3 many times, then plot correlation.


setmydirs_detectv3;
% show ppant numbers:
cd(procdatadir)
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%   


%% concat data:
%preallocate storage:

pidx1=ceil(linspace(1,100,21)); % length n-1
pidx2=ceil(linspace(1,100,41));% 
useindexes = {pidx1, pidx2};
gaittypes = {'single gait' , 'double gait'};
%
pidx=pidx2;
 mdiff = round(mean(diff(pidx)./2));
    xvec = pidx(1:end-1) + mdiff;
    %%
for ippant =1%:nsubs
    cd(procdatadir)    
    load(pfols(ippant).name, 'subjID', 'trial_summaryTable',...
        'gait_ts_gData', 'gait_ts_resamp','gait_ts_eyedirY_resamp', ...
        'gait_ts_eyedirZ_resamp', 'gait_ts_eyeposY_resamp', 'gait_ts_eyeposZ_resamp',...
        'doubgait_ts_resamp');
  
    
    disp(['performing split-half on subject ' num2str(ippant)]);
    
    
    %% also calculate binned versions
    % These take the mean over an index range, per gait cycle position
    
    nGait = 2;  % use to have option for 1 or 2.
    useL =~isnan(trial_summaryTable.trgO_gPcnt_LRL);
        useR =~isnan(trial_summaryTable.trgO_gPcnt_RLR);
       
        pidx=useindexes{nGait}; % indices for binning (single or double gc).
        ppantData= trial_summaryTable;
    

% for different permutations, we will subselect from the data:
% Divide dataset into two halves

            %% %%%%%%%%%%%%%%%%
            % Step through different data types :
            
            % note that now, we will sample randomly from all
            % positions, based on the number in real data.
            
            % 3 main analysis (for behaviour).
            %- RT relative to target gait position
            %- RT relative to response gait position
            %- Acc relative to target gait position.
            % nb: acc relative to response, is moot, since incorrect responses
            %(False Alarms) were so few.
            %%%%%%%%%%%%%%%%%%
            % for each of these three types, define the search values
            % to shuffle over:
            
           
            searchPosIDX = {[ppantData.trgO_gPcnt_LRL(useL); ppantData.trgO_gPcnt_RLR(useR)], ... % trg onset as %               
                [ppantData.trgO_gPcnt_LRL(useL);ppantData.trgO_gPcnt_RLR(useR)]}; % targ onset as %

           
             %and define the DV of interest:
            searchPosDVs = {[ppantData.clickRT(useL);ppantData.clickRT(useR)], ... %RT, RT, Acc.
                [ppantData.targCor(useL); ppantData.targCor(useR)]};

            dv=nan;
            idx=nan;
%             accTable=table(dv, idx);
%             rtTable=table();
%             respTable=table();

            for itype=1:2

                collectDVs=[]
                collectIdx=[];
                %preallocate:
                outgoing = nan(1,length(pidx)); %

                outgoingALL=[];
                allpos = searchPosIDX{itype};
                allDVs = searchPosDVs{itype};
                
                
                % for actual gait pos, select reldata
                for ip=1:pidx(end)
                    useb = find(allpos==ip);
                    outgoingALL(ip).d= allDVs(useb);
                    
                    collectDVs= [collectDVs; allDVs(useb)];
                    collectIdx= [collectIdx; repmat(ip, [length(useb),1])]; 
                end
                
                switch itype
                    case 1
                      
                        targOns_RTs_ALL = table(collectDVs, collectIdx);
                   
                    case 2
                        targOns_Acc_ALL= table(collectDVs, collectIdx);
                end
            end % itype
            


            %% %%%%%%%%%%%%%%%%%%%
            % Step through all data types, performing binning based on
            % length of pidx.
            %%%%%%%%%%%%%%%%%%%
            %trg onset counts
            %resp onset counts
            %acc per trg onset
            %acc per resp onset
            %rt per trg onset
            %rt per resp onset
            

% bin per datatype, based on split.


            usedata= {trgOtoAvg, respOtoAvg, ... % using ver with NaNs removed
                targOns_Acc,targOns_RTs, respOns_RTs};
            
            for itype= 1:5
                
                %% for all types, take mean over bin indices.
                out_bin=nan(1,length(pidx)-1);
                datatobin = usedata{itype};
                for ibin=1:length(pidx)-1
                    idx = pidx(ibin):pidx(ibin+1);
                    
                    out_bin(ibin) = mean(datatobin(idx), 'omitnan');
                    
                end
                
                switch itype
                    case 1
                        trgCounts_bin = out_bin;
                    case 2
                        respCounts_bin = out_bin;
                    case 3
                        targOns_Acc_bin = out_bin;                   
                        
                    case 4
                        targOns_RT_bin= out_bin;
                    case 5
                        respOns_RT_bin = out_bin;
                end
                
            end

%% now we can perform the subj fit:
for irep= 1

n = size(ppantData, 1);
idx = randperm(n);
half1 = idx(1:n/2);
half2 = idx(n/2+1:end);

usesubstack= {half1, half2};
for isplit= 1:2

    uset=usesubstack{isplit};


    Hzspace = [0.01:.2:10];
    
    %best fit (use model params below):
    gM= targOns_Acc_bin;
    f = fit(xvec', gM', 'fourier1'); %unbounded
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

testfreqs=[2,4];
   for ifreq=[1,2]; 
   
        % include period and Rsquared
        
        testw = 2*pi*testfreqs(ifreq)/xvec(end);
        
        CoeffBounds.w(1) = testw;
        CoeffBounds.w(2) = testw;
        
        %update fit opts settings
        
        %Update Fit Options setting.
        FitOpts = fitoptions('Method','NonlinearLeastSquares','Lower',table2array(CoeffBounds(1,:)),...
            'Upper',table2array(CoeffBounds(2,:)));
        
        %first test this period on observed data
        %set last coefficient value in the fourier model (w) to this value:
        
        [tmpf,gof] = fit(xvec', gM', 'fourier1', FitOpts);
        
        % how good/bad was the fit?
            allFits(irep,isplit,ifreq) = gof.rsquare;

            %% %%%%%%%%%%%%%%%%%%%%%%%%% now save:
            
   end   %ifreq
        end % isplit
        disp(['finished rep ' num2str(irep)])
    end % irep.
end % ippant
