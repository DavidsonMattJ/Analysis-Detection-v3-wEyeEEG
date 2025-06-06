% j5A_createNullGFX_sacc
%% similar to the concat job (j4), except we disregard gait position information, 
 % and instead sample uniformly from all positions to create a null distribution.
%  [nbins x nPerm] per subject and datatype.



setmydirs_detectv3;
%% show ppant numbers:
cd(procdatadir);
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%


%% concat data:
%preallocate storage:

gaittypes = {'single gait' , 'double gait'};
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %preallocate storage:
   
     %use same spacing as observed data, for the bins:
     %add extra fields for the shuffled data.
     cd([procdatadir filesep 'GFX'])
    load('GFX_Data_inGaits', 'GFX_SaccadePosData','pidx1', 'pidx2');
    
%     GFX_TargPos_nullData=[]; % store the results per subj, for stats.
%     GFX_RespPos_nullData=[];
    
    nPerm = 1000; % how many times to perform the resampling?
    
    for ippant =1:length(pfols)
        
        cd(procdatadir)    %%load data from import job.
        load(pfols(ippant).name, ...
            'trial_summaryTable','subjID');
        ppantData= saccadeSummaryTable;

%       
       
        disp(['constructing null distributions for ' subjID]);
 
        %%

        useindexes= {pidx1, pidx2};
        for nGait=2%1:2  % use to have option for 1 or 2.

            pidx=useindexes{nGait}; % indices for binning (single or double gc).

            if nGait==1 % single step.

                % first retrieve index for separate gaits (L/Right feet).
                useL= strcmp(ppantData.saccOns_gFoot, 'LR');
                useR= strcmp(ppantData.saccOns_gFoot, 'RL');
                useAll = 1:length(useL);

                % split by foot (LR, RL, all combined)
            else

                % As above, just use all trials per speed.
                useall= 1:length(useL);


            end
            trialstoIndex ={find(useL), find(useR), useAll};



        for iLR=3% 1:3 % Left Right, both
           
               uset=find(trialstoIndex{iLR});
                
                

               if nGait==1
                   searchPosIDX = ppantData.saccOnset_gPcnt(uset); % targ onset as %

               else

                   searchPosIDX = [ppantData.saccOnset_gPcnt_step1inStride(uset); ppantData.saccOnset_gPcnt_step2inStride(uset)];

               end


               useposdata = searchPosIDX;

               % raw hit counts, per gaitsample:
               real_counts= histcounts(useposdata, pidx(end));

               %critical step -convert 0 to nans. this avoids averaging
               %problems (as above, in data anaysis)

               real_counts(real_counts==0)=NaN;
               %now create a null, by shuffling these counts per perm.

               % preallocate:
               outgoingnull = nan(nPerm, pidx(end));
               for iperm = 1:nPerm


                   a = rand(1,length(real_counts));
                   [~, i] = sort(a);
                   outgoingnull(iperm,:) =real_counts(i);
               end


               saccOns_Counts_null = outgoingnull;


                %% %%%%%%%%%%%%%%%%%%%
                % Step through all shuffled data, storing binned versions, per
                % permutation previously performed.
                %%%%%%%%%%%%%%%%%%%
               
                binmedata = {saccOns_Counts_null};%, respOns_RTs_null, targOns_Acc_null,...
%                     targOns_Counts_null, respOns_Counts_null};
                
                for itype = 1%:length(binmedata)
                   
                   usedata = binmedata{itype};                   
                   outgoingbinned=nan(nPerm, length(pidx)-1);
                   
                   for ibin=1:length(pidx)-1
                       idx = pidx(ibin):pidx(ibin+1);
                       
                       outgoingbinned(:,ibin) = nanmean(usedata(:,idx),2);
                       
                   end
                    
                   %rename this time. to avoid headaches
                   switch itype
                       case 1
                           saccCounts_bin_shuff= outgoingbinned;
%                        case 2
%                            respRT_bin_shuff= outgoingbinned;
%                        case 3
%                            trgAcc_bin_shuff= outgoingbinned;
%                        case 4
%                            trgCounts_bin_shuff = outgoingbinned;
%                        case 5
%                            respCounts_bin_shuff = outgoingbinned;
                   end
                end % per type
                
               %% store for all ppants, preparing to save:
                gaitnames = {'gc_', 'doubgc_'}; % change fieldname based on nGaits:
                
                %add binned, shuffled vers:
                GFX_SaccadePosData(ippant,iLR).([gaitnames{nGait} 'sacc_all_binned_counts_shuff']) = saccCounts_bin_shuff;
                

                
            end % iLR
        end
 

    disp(['finished creation of null sacc onset data for ppant ' num2str(ippant)]);
    end % ppant
    
    %%
    cd([procdatadir filesep 'GFX']);
    
    disp('resaving GFX data with shuffle attached...');
    save('GFX_Data_inGaits', ...
      'GFX_SaccadePosData','-append');